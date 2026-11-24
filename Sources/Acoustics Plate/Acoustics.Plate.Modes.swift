public import Acoustics
public import Eigen
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
    @usableFromInline
    internal struct System {
        @usableFromInline let stiffness: Sparse.Symmetric<Double>
        @usableFromInline let mass: Sparse.Symmetric<Double>
        /// For each vertex, the unknown its deflection is, or nil where it is held.
        @usableFromInline let deflection: [Int?]
        /// Each face's corners, counterclockwise, as vertex indices, and the element over them.
        @usableFromInline let faces: [[Int]]
        @usableFromInline let elements: [Element]

        @usableFromInline var order: Int { stiffness.order }
    }
}

extension Acoustics.Plate.System {
    @usableFromInline
    init(_ plate: Acoustics.Plate) throws(Acoustics.Plate.Error) {
        let points = plate.mesh.vertices
        var faces: [[Int]] = []
        var elements: [Acoustics.Plate.Element] = []
        for (index, face) in plate.mesh.faces.enumerated() {
            var corners = face.vertices.map { Int(bitPattern: $0) }
            guard corners.count == 4 else { throw .face(index) }
            var element = Acoustics.Plate.Element(x: corners.map { points[$0].x }, y: corners.map { points[$0].y })
            if element.determinant < 0 {
                corners.reverse()
                element = Acoustics.Plate.Element(x: corners.map { points[$0].x }, y: corners.map { points[$0].y })
            }
            guard element.determinant > 0 else { throw .face(index) }
            faces.append(corners)
            elements.append(element)
        }
        let held = Self.held(plate.edge, exceptions: plate.exceptions, points: points, faces: faces)

        // Number the unknowns: deflection, then the two rotations in each vertex's frame.
        var unknown = [Int](repeating: -1, count: 3 * points.count)
        var order = 0
        for vertex in points.indices {
            for component in 0..<3 where !held.fixed[vertex].contains(component) {
                unknown[3 * vertex + component] = order
                order += 1
            }
        }

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

        var stiffness: [Sparse.Symmetric<Double>.Entry] = []
        var mass: [Sparse.Symmetric<Double>.Entry] = []
        stiffness.reserveCapacity(78 * faces.count)
        mass.reserveCapacity(78 * faces.count)
        for (corners, element) in zip(faces, elements) {
            var (k, m) = element.matrices(thickness: h, density: material.density.underlying, q: q, s: s)
            for (local, vertex) in corners.enumerated() {
                guard let tangent = held.tangent[vertex] else { continue }
                Self.rotate(&k, corner: local, tangent: tangent)
                Self.rotate(&m, corner: local, tangent: tangent)
            }
            var map = [Int](repeating: -1, count: 12)
            for (local, vertex) in corners.enumerated() {
                for component in 0..<3 { map[3 * local + component] = unknown[3 * vertex + component] }
            }
            for a in 0..<12 where map[a] >= 0 {
                for b in a..<12 where map[b] >= 0 {
                    let row = Ordinal(UInt(map[a]))
                    let column = Ordinal(UInt(map[b]))
                    stiffness.append(Sparse.Symmetric.Entry(row, column, k[12 * a + b]))
                    mass.append(Sparse.Symmetric.Entry(row, column, m[12 * a + b]))
                }
            }
        }
        do throws(Sparse.Error) {
            self.stiffness = try Sparse.Symmetric(order: Cardinal(UInt(order)), entries: stiffness)
            self.mass = try Sparse.Symmetric(order: Cardinal(UInt(order)), entries: mass)
        } catch {
            preconditionFailure("every unknown is numbered below the order")
        }
        self.deflection = points.indices.map { unknown[3 * $0] >= 0 ? unknown[3 * $0] : nil }
        self.faces = faces
        self.elements = elements
    }

    /// Turns a corner's two rotations from (βx, βy) into (along the tangent, along its normal): with T the rotation,
    /// the matrix becomes Tᵀ A T.
    @usableFromInline
    static func rotate(_ matrix: inout [Double], corner: Int, tangent: (x: Double, y: Double)) {
        let first = 3 * corner + 1
        let second = first + 1
        // β = ψt t + ψn n, with n = (−ty, tx): column first of T is t, column second is n.
        for row in 0..<12 {
            let a = matrix[12 * row + first]
            let b = matrix[12 * row + second]
            matrix[12 * row + first] = a * tangent.x + b * tangent.y
            matrix[12 * row + second] = -a * tangent.y + b * tangent.x
        }
        for column in 0..<12 {
            let a = matrix[12 * first + column]
            let b = matrix[12 * second + column]
            matrix[12 * first + column] = a * tangent.x + b * tangent.y
            matrix[12 * second + column] = -a * tangent.y + b * tangent.x
        }
    }

    /// Which of each vertex's unknowns the edge and its exceptions hold, and the tangent of the outline where only
    /// the rotation about the edge's normal is held.
    ///
    /// The outline is the mesh's edges that only one face uses. For `.supported`, a vertex where the outline turns
    /// by no more than 45° keeps its rotation about the outline, along the tangent that bisects its two edges; any
    /// other vertex of the outline, such as a corner, has both rotations held.
    @usableFromInline
    static func held(
        _ edge: Acoustics.Plate.Edge,
        exceptions: [Int: Acoustics.Plate.Edge],
        points: [Point<2, Double>],
        faces: [[Int]]
    ) -> (fixed: [Set<Int>], tangent: [(x: Double, y: Double)?]) {
        var fixed = [Set<Int>](repeating: [], count: points.count)
        var tangent = [(x: Double, y: Double)?](repeating: nil, count: points.count)
        guard edge != .free || exceptions.values.contains(where: { $0 != .free }) else { return (fixed, tangent) }
        var uses: [Int: Int] = [:]
        func key(_ a: Int, _ b: Int) -> Int { Swift.min(a, b) * points.count + Swift.max(a, b) }
        for corners in faces {
            for i in 0..<corners.count {
                uses[key(corners[i], corners[(i + 1) % corners.count]), default: 0] += 1
            }
        }
        var directions = [[(x: Double, y: Double)]](repeating: [], count: points.count)
        for corners in faces {
            for i in 0..<corners.count {
                let a = corners[i]
                let b = corners[(i + 1) % corners.count]
                guard uses[key(a, b)] == 1 else { continue }
                let dx = points[b].x - points[a].x
                let dy = points[b].y - points[a].y
                let length = (dx * dx + dy * dy).squareRoot()
                directions[a].append((dx / length, dy / length))
                directions[b].append((dx / length, dy / length))
            }
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
                    let length = (tx * tx + ty * ty).squareRoot()
                    tangent[vertex] = (tx / length, ty / length)
                    fixed[vertex] = [0, 1]
                } else {
                    fixed[vertex] = [0, 1, 2]
                }
            }
        }
        return (fixed, tangent)
    }
}

extension Acoustics.Plate.System {
    /// The mode of an eigenpair: its frequency, its shape scaled to a largest deflection of one, its mass at that
    /// scale, and its radiation into `air`.
    @usableFromInline
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
    @usableFromInline
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
