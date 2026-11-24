import Bessel
import Numeric
import Testing

@testable import Acoustics_Plate

extension `Acoustics.Plate Tests` {
    /// Radiation into a half-space through an infinite rigid baffle, against closed forms.
    @Suite struct Radiation {}
}

extension `Acoustics.Plate Tests`.Radiation {
    @Test(arguments: [0.25, 0.5, 1, 2, 4])
    func `A disk moving as one radiates as a baffled piston, 1 − J₁(2ka) / ka`(ka: Double) throws {
        let radius = 0.1
        let plate = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.disk(radius: radius, cells: 12),
            thickness: Space.Physical.metres(0.001),
            material: `Acoustics.Plate Tests`.aluminium,
            edge: .free
        )
        let system = try Acoustics.Plate.System(plate)
        let shape = [Double](repeating: 1, count: plate.mesh.vertices.count)
        let found = system.radiation(shape, wavenumber: ka / radius)
        let exact = 1 - 2 * Bessel.First<Double>(order: 1).value(at: 2 * ka) / (2 * ka)
        #expect(
            abs(found.efficiency - exact) < 0.01 * Swift.max(exact, 0.1), "found \(found.efficiency), piston \(exact)")
        // The mesh's outline is a polygon inside the circle.
        #expect(abs(found.area - Double.pi * radius * radius) < 0.005 * Double.pi * radius * radius)
    }

    @Test
    func `A simply supported panel's fundamental radiates at low frequency as Wallace's 32 k²ab / π⁵`() throws {
        let a = 0.4
        let b = 0.25
        let plate = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.rectangle(a, b, columns: 32, rows: 20),
            thickness: Space.Physical.metres(0.003),
            material: Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter),
            edge: .supported
        )
        let system = try Acoustics.Plate.System(plate)
        let shape = plate.mesh.vertices.map { Double.sin(.pi * $0.x / a) * Double.sin(.pi * $0.y / b) }
        let k = 0.05 / a
        // Wallace (1972), J. Acoust. Soc. Am. 51, 946: σ ≈ 32 k² a b / (π⁵ m² n²) for odd m, n when ka ≪ 1.
        let wallace = 32 * k * k * a * b / (Double.pi * .pi * .pi * .pi * .pi)
        let found = system.radiation(shape, wavenumber: k)
        #expect(abs(found.efficiency - wallace) / wallace < 0.01, "found \(found.efficiency), Wallace \(wallace)")
        // ∫∫ sin(πx/a) sin(πy/b) = 4ab / π².
        #expect(abs(found.area - 4 * a * b / (.pi * .pi)) / (4 * a * b / (.pi * .pi)) < 0.005)
    }

    @Test
    func `An even mode pushes out no net volume and radiates far less than an odd one`() throws {
        let plate = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.rectangle(0.4, 0.25, columns: 32, rows: 20),
            thickness: Space.Physical.metres(0.003),
            material: Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter),
            edge: .supported
        )
        let modes = try Acoustics.Plate.Modes(plate, count: 2)
        #expect(abs(modes[1].area.underlying) < 1e-9)
        #expect(modes[1].efficiency < modes[0].efficiency)
        #expect(modes[0].area.underlying > 0)
    }
}
