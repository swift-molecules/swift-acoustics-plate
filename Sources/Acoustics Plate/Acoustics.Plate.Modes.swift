public import Acoustics
public import Eigen
internal import Galerkin
internal import Numeric

extension Acoustics.Plate {
    /// A plate's lowest modes, computed once from its shape and material and kept as a value, with each mode's
    /// radiation into the air given.
    public struct Modes {
        /// The plate they belong to.
        public let plate: Acoustics.Plate
        /// The air they radiate into.
        public let air: Acoustics.Air
        /// The modes, lowest first; a repeated frequency appears once for each mode that has it.
        public let elements: [Mode]

        /// The `count` lowest modes of `plate`, radiating into `air`.
        ///
        /// A free plate's first three modes are its rigid motions, at zero frequency to rounding.
        public init(_ plate: Acoustics.Plate, count: Int, air: Acoustics.Air = .room) throws(Acoustics.Plate.Error) {
            let system = try Acoustics.Plate.System(plate)
            guard count <= system.order else { throw .count }
            let pairs: [Eigen.Pair<Double>]
            do throws(Eigen.Error) {
                let pencil = try Eigen.Pencil(system.stiffness, system.mass)
                pairs = try pencil.lowest(Cardinal(UInt(count)))
            } catch {
                throw .eigen(error)
            }
            self.plate = plate
            self.air = air
            self.elements = pairs.map { system.mode($0, air: air) }
        }
    }
}

extension Acoustics.Plate.Modes: Hashable {}
extension Acoustics.Plate.Modes: Sendable {}

extension Acoustics.Plate.Modes: RandomAccessCollection {
    public var startIndex: Int { elements.startIndex }
    public var endIndex: Int { elements.endIndex }
    public subscript(position: Int) -> Acoustics.Plate.Mode { elements[position] }
}

extension Acoustics.Plate {
    /// A plate's stiffness and mass over its unknowns that are not held, and what turns an eigenvector back into a
    /// mode.
    internal struct System {
        let stiffness: Sparse.Symmetric<Double>
        let mass: Sparse.Symmetric<Double>
        /// For each vertex, the unknown its deflection is, or nil where it is held.
        let deflection: [Int?]
        /// Each face's corners, counterclockwise, as vertex indices, and the element over them.
        let faces: [[Int]]
        let elements: [Element]

        var order: Int { stiffness.order }
    }
}

extension Acoustics.Plate.System {
    init(_ plate: Acoustics.Plate) throws(Acoustics.Plate.Error) {
        let points = plate.mesh.vertices
        var faces: [[Int]] = []
        var elements: [Acoustics.Plate.Element] = []
        for (index, face) in plate.mesh.faces.enumerated() {
            var corners = face.vertices.map { Int(bitPattern: $0) }
            guard corners.count == 4 else { throw .face(index) }
            var element = Acoustics.Plate.Element(corners.map { points[$0] })
            if element.determinant < 0 {
                corners.reverse()
                element = Acoustics.Plate.Element(corners.map { points[$0] })
            }
            guard element.determinant > 0 else { throw .face(index) }
            faces.append(corners)
            elements.append(element)
        }
        let held = Self.held(plate.edge, exceptions: plate.exceptions, points: points, outline: plate.mesh.outline)

        // The unknowns: deflection, then the two rotations in each vertex's frame.
        let unknowns = Galerkin.Unknowns(vertices: points.count, components: 3, held: held.fixed)

        let material = plate.material
        let ex = material.young.x.underlying
        let ey = material.young.y.underlying
        let nxy = material.poisson.underlying
        let nyx = nxy * ey / ex
        let denominator = 1 - nxy * nyx
        let q = (
            q11: ex / denominator, q12: nxy * ey / denominator, q22: ey / denominator, q66: material.shear.xy.underlying
        )
        let h = plate.thickness.underlying
        let s = (xz: 5.0 / 6 * h * material.shear.xz.underlying, yz: 5.0 / 6 * h * material.shear.yz.underlying)

        var stiffness = Galerkin.Assembly(unknowns, rotations: held.rotations)
        var mass = stiffness
        for (corners, element) in zip(faces, elements) {
            let (k, m) = element.matrices(thickness: h, density: material.density.underlying, q: q, s: s)
            stiffness.add(k, at: corners)
            mass.add(m, at: corners)
        }
        self.stiffness = stiffness.matrix
        self.mass = mass.matrix
        self.deflection = points.indices.map { unknowns[$0, 0] }
        self.faces = faces
        self.elements = elements
    }

    /// Which of each vertex's unknowns the edge and its exceptions hold, and where only the rotation about the edge's
    /// normal is held, the rotations turned into the outline's tangent and normal.
    ///
    /// The outline is the mesh's `outline`, the edges that only one face uses. For `.supported`, a vertex where the
    /// outline turns by no more than 45° keeps its rotation about the outline, along the tangent that bisects its two
    /// edges; any other vertex of the outline, such as a corner, has both rotations held.
    static func held(
        _ edge: Acoustics.Plate.Edge,
        exceptions: [Int: Acoustics.Plate.Edge],
        points: [Point<2, Double>],
        outline: [Segment<Index<Point<2, Double>>>]
    ) -> (fixed: [Int: Set<Int>], rotations: [Int: Galerkin.Rotation<Double>]) {
        var fixed: [Int: Set<Int>] = [:]
        var rotations: [Int: Galerkin.Rotation<Double>] = [:]
        guard edge != .free || exceptions.values.contains(where: { $0 != .free }) else { return (fixed, rotations) }
        var directions = [[(x: Double, y: Double)]](repeating: [], count: points.count)
        for segment in outline {
            let a = Int(bitPattern: segment.start)
            let b = Int(bitPattern: segment.end)
            let dx = points[b].x - points[a].x
            let dy = points[b].y - points[a].y
            let length = (dx * dx + dy * dy).squareRoot()
            directions[a].append((dx / length, dy / length))
            directions[b].append((dx / length, dy / length))
        }
        let straight = 0.5.squareRoot()
        for vertex in points.indices where !directions[vertex].isEmpty {
            switch exceptions[vertex] ?? edge {
            case .free:
                break
            case .pinned:
                fixed[vertex] = [0]
            case .clamped:
                fixed[vertex] = [0, 1, 2]
            case .supported:
                let d = directions[vertex]
                let turn = d.count == 2 ? d[0].x * d[1].x + d[0].y * d[1].y : 0
                if d.count == 2 && Swift.abs(turn) >= straight {
                    let sign: Double = turn < 0 ? -1 : 1
                    let tx = d[0].x + sign * d[1].x
                    let ty = d[0].y + sign * d[1].y
                    // (βx, βy) = ψt t + ψn n: ψt, the first after the turn, is held with w.
                    rotations[vertex] = Galerkin.Rotation(1, 2, along: (x: tx, y: ty))
                    fixed[vertex] = [0, 1]
                } else {
                    fixed[vertex] = [0, 1, 2]
                }
            }
        }
        return (fixed, rotations)
    }
}

extension Acoustics.Plate.System {
    /// The mode of an eigenpair: its frequency, its shape scaled to a largest deflection of one, its mass at that
    /// scale, and its radiation into `air`.
    func mode(_ pair: Eigen.Pair<Double>, air: Acoustics.Air) -> Acoustics.Plate.Mode {
        var shape = deflection.map { $0.map { pair.vector[$0] } ?? 0 }
        var peak = 0.0
        for value in shape where Swift.abs(value) > Swift.abs(peak) { peak = value }
        if peak != 0 {
            for index in shape.indices { shape[index] /= peak }
        }
        let frequency = Swift.max(pair.value, 0).squareRoot() / (2 * Double.pi)
        // The vector is B-normalized, so at a largest deflection of one the modal mass is 1 / peak².
        let mass = peak != 0 ? 1 / (peak * peak) : 0
        let radiation = radiation(shape, wavenumber: 2 * Double.pi * frequency / air.speed)
        return Acoustics.Plate.Mode(
            frequency: Frequency.hertz(frequency),
            shape: shape,
            mass: Mass.kilograms(mass),
            area: Area<Space.Physical>.amount(radiation.area),
            efficiency: radiation.efficiency
        )
    }

    /// The signed swept area and the baffled radiation efficiency of a shape vibrating at wavenumber `k`.
    ///
    /// Each face is a point source at its centroid of strength ∫ w dA (the Rayleigh integral, one point per face),
    /// so the radiated power over ρ₀ c is (k² / 4π) Σᵢ Σⱼ qᵢ qⱼ sin(k rᵢⱼ) / (k rᵢⱼ) per square of the angular
    /// frequency, and the efficiency is that over ½ ∫ w² dA.
    func radiation(_ shape: [Double], wavenumber k: Double) -> (area: Double, efficiency: Double) {
        var strength = [Double](repeating: 0, count: faces.count)
        var cx = strength
        var cy = strength
        var square = 0.0
        for (index, (corners, element)) in zip(faces, elements).enumerated() {
            let moments = element.moments
            var q = 0.0
            for i in 0..<4 { q += moments.weights[i] * shape[corners[i]] }
            strength[index] = q
            cx[index] = moments.x / moments.area
            cy[index] = moments.y / moments.area
            let overlap = element.overlap
            for i in 0..<4 {
                for l in 0..<4 { square += overlap[4 * i + l] * shape[corners[i]] * shape[corners[l]] }
            }
        }
        let area = strength.reduce(0, +)
        guard k > 0, square > 0 else { return (area, 0) }
        var sum = 0.0
        for i in strength.indices {
            sum += strength[i] * strength[i]
            var j = i + 1
            while j < strength.count {
                let dx = cx[i] - cx[j]
                let dy = cy[i] - cy[j]
                let kr = k * (dx * dx + dy * dy).squareRoot()
                sum += 2 * strength[i] * strength[j] * (kr > 0 ? Double.sin(kr) / kr : 1)
                j += 1
            }
        }
        return (area, k * k / (4 * Double.pi) * sum / (0.5 * square))
    }
}
