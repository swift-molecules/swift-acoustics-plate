@testable import Acoustics_Plate
import Testing

extension `Acoustics.Plate Tests` {
    /// What holds for every plate.
    @Suite struct Laws {}
}

extension `Acoustics.Plate Tests`.Laws {
    static let panel = Acoustics.Plate(
        mesh: `Acoustics.Plate Tests`.rectangle(0.4, 0.25, columns: 16, rows: 10),
        thickness: Space.Physical.metres(0.003),
        material: Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter),
        edge: .clamped
    )

    @Test
    func `Modes are ascending, their shapes peak at one, held vertices do not move, and the mass is positive`() throws {
        let modes = try Acoustics.Plate.Modes(Self.panel, count: 8)
        #expect(modes.count == 8)
        for (lower, higher) in zip(modes, modes.dropFirst()) {
            #expect(lower.frequency.underlying <= higher.frequency.underlying)
        }
        for mode in modes {
            #expect(mode.shape.max() == 1)
            #expect(mode.shape.min()! >= -1)
            #expect(mode.mass.underlying > 0)
            for (point, value) in zip(Self.panel.mesh.vertices, mode.shape)
            where point.x == 0 || point.y == 0 || point.x == 0.4 || point.y == 0.25 {
                #expect(value == 0)
            }
        }
    }

    @Test
    func `The modes are a value: computed twice they are equal`() throws {
        #expect(try Acoustics.Plate.Modes(Self.panel, count: 4) == Acoustics.Plate.Modes(Self.panel, count: 4))
    }

    @Test
    func `A plate scaled by two in every length rings at half the frequency`() throws {
        let doubled = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.rectangle(0.8, 0.5, columns: 16, rows: 10),
            thickness: Space.Physical.metres(0.006),
            material: Self.panel.material,
            edge: .clamped
        )
        let small = try Acoustics.Plate.Modes(Self.panel, count: 4)
        let large = try Acoustics.Plate.Modes(doubled, count: 4)
        for (s, l) in zip(small, large) {
            #expect(abs(2 * l.frequency.underlying - s.frequency.underlying) < 1e-9 * s.frequency.underlying)
            #expect(abs(l.mass.underlying - 8 * s.mass.underlying) < 1e-9 * l.mass.underlying)
        }
    }

    @Test
    func `Faces listed clockwise give the same modes as counterclockwise`() throws {
        let reversed = Acoustics.Plate(
            mesh: Self.panel.mesh.reversed,
            thickness: Self.panel.thickness,
            material: Self.panel.material,
            edge: Self.panel.edge
        )
        let a = try Acoustics.Plate.Modes(Self.panel, count: 3)
        let b = try Acoustics.Plate.Modes(reversed, count: 3)
        for (x, y) in zip(a, b) {
            #expect(abs(x.frequency.underlying - y.frequency.underlying) < 1e-9 * x.frequency.underlying)
        }
    }

    @Test
    func `Clamping, then supporting, then pinning, then freeing the edge lowers every frequency`() throws {
        var previous: [Double]? = nil
        for edge in [Acoustics.Plate.Edge.clamped, .supported, .pinned, .free] {
            var plate = Self.panel
            plate.edge = edge
            let found = `Acoustics.Plate Tests`.frequencies(try Acoustics.Plate.Modes(plate, count: 4))
            if let previous {
                for (p, f) in zip(previous, found) { #expect(f <= p * (1 + 1e-9)) }
            }
            previous = found
        }
    }

    @Test
    func `Asking for more modes than the plate has fails`() {
        let tiny = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.rectangle(0.1, 0.1, columns: 1, rows: 1),
            thickness: Space.Physical.metres(0.001),
            material: `Acoustics.Plate Tests`.aluminium,
            edge: .free
        )
        #expect(throws: Acoustics.Plate.Error.count) { try Acoustics.Plate.Modes(tiny, count: 13) }
    }

    @Test
    func `A triangle is refused`() throws {
        let mesh = try Mesh(
            vertices: [Point<2, Double>(x: 0, y: 0), Point(x: 1, y: 0), Point(x: 0, y: 1)],
            faces: [Polygon(0, 1, 2)]
        )
        let plate = Acoustics.Plate(
            mesh: mesh,
            thickness: Space.Physical.metres(0.001),
            material: `Acoustics.Plate Tests`.aluminium,
            edge: .free
        )
        #expect(throws: Acoustics.Plate.Error.face(0)) { try Acoustics.Plate.Modes(plate, count: 1) }
    }
}

extension `Acoustics.Plate Tests` {
    /// The woods are the USDA Wood Handbook's (2021) tables, as the source states them.
    @Suite struct Wood {}
}

extension `Acoustics.Plate Tests`.Wood {
    @Test
    func `Sitka spruce: E_L is Table 5–3a's 10,800 MPa raised 10%, the rest Table 5–1's ratios, density from G₁₂ = 0.40`() {
        let wood = Acoustics.Plate.Material.Wood.sitkaSpruce
        let e = 11_880e6
        #expect(abs(wood.young.longitudinal.underlying - e) < 1)
        #expect(abs(wood.young.radial.underlying - 0.078 * e) < 1)
        #expect(abs(wood.young.tangential.underlying - 0.043 * e) < 1)
        #expect(abs(wood.shear.longitudinalRadial.underlying - 0.064 * e) < 1)
        #expect(abs(wood.shear.longitudinalTangential.underlying - 0.061 * e) < 1)
        #expect(abs(wood.shear.radialTangential.underlying - 0.003 * e) < 1)
        #expect(wood.poisson.longitudinalRadial.underlying == 0.372)
        #expect(wood.poisson.longitudinalTangential.underlying == 0.467)
        #expect(abs(wood.density.underlying - 448) < 1e-9)
    }

    @Test
    func `Quartersawn puts the radial axis across the plate, flatsawn the tangential`() {
        let wood = Acoustics.Plate.Material.Wood.engelmannSpruce
        let quarter = Acoustics.Plate.Material(wood, cut: .quarter)
        let flat = Acoustics.Plate.Material(wood, cut: .flat)
        #expect(quarter.young.y == wood.young.radial)
        #expect(flat.young.y == wood.young.tangential)
        #expect(quarter.shear.xy == wood.shear.longitudinalRadial)
        #expect(flat.shear.xy == wood.shear.longitudinalTangential)
        #expect(quarter.shear.yz == wood.shear.radialTangential)
        #expect(flat.shear.yz == wood.shear.radialTangential)
        #expect(quarter.poisson == wood.poisson.longitudinalRadial)
        #expect(flat.poisson == wood.poisson.longitudinalTangential)
    }
}
