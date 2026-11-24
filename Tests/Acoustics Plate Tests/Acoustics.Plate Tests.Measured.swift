import Testing

@testable import Acoustics_Plate

extension `Acoustics.Plate Tests` {
    /// Spruce boards whose frequencies were measured, computed from the constants identified for them.
    ///
    /// Ducceschi, Duran, Tahvanainen and Ausiello (2024), "A method to estimate the rectangular orthotropic plate
    /// elastic constants using least-squares and Chladni patterns", Applied Acoustics 220, 109949 (open access,
    /// CC BY; arXiv:2311.06935): the boards' sizes and densities (Table 1 and section 4.2), their measured
    /// frequencies (Tables 8 and 9), and the mean E_x, E_y and G_xy identified from them (Table 10). ν_xy and the
    /// transverse shear moduli G_xz and G_yz are the ones their Table 1 assumes. Their own finite-element model
    /// (COMSOL) with the same constants is within 2% of the measurements (Tables 11 and 12).
    @Suite struct Measured {}
}

extension `Acoustics.Plate Tests`.Measured {
    static func spruce(density: Double, x: Double, y: Double, shear: Double) -> Acoustics.Plate.Material {
        Acoustics.Plate.Material(
            density: Density.amount(density),
            young: .init(x: Modulus.Young.pascals(x), y: Modulus.Young.pascals(y)),
            shear: .init(
                xy: Modulus.Shear.pascals(shear), xz: Modulus.Shear.pascals(0.62e9), yz: Modulus.Shear.pascals(0.023e9)),
            poisson: Poisson.amount(0.51)
        )
    }

    @Test
    func `A cantilever Finnish spruce board rings within 3% of its six measured frequencies`() throws {
        let a = 0.223
        let b = 0.114
        let mesh = `Acoustics.Plate Tests`.rectangle(a, b, columns: 40, rows: 20)
        var clamped: [Int: Acoustics.Plate.Edge] = [:]
        for (index, point) in mesh.vertices.enumerated() where point.x == 0 { clamped[index] = .clamped }
        let plate = Acoustics.Plate(
            mesh: mesh,
            thickness: Space.Physical.metres(0.003),
            material: Self.spruce(density: 473.9, x: 13.0e9, y: 702e6, shear: 668e6),
            edge: .free,
            exceptions: clamped
        )
        let found = `Acoustics.Plate Tests`.frequencies(try Acoustics.Plate.Modes(plate, count: 6))
        // Table 8: (1,0), (1,1), (2,0), (1,2), (2,1), (2,2).
        let measured = [52.0, 98, 310, 337, 398, 637]
        for (f, m) in zip(found, measured) {
            #expect(abs(f - m) / m < 0.03, "computed \(f) Hz, measured \(m) Hz")
        }
    }

    @Test
    func `A clamped Fiemme Valley red spruce guitar top rings within 3% of its six measured frequencies`() throws {
        let plate = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.rectangle(0.60, 0.412, columns: 48, rows: 33),
            thickness: Space.Physical.metres(0.003),
            material: Self.spruce(density: 399, x: 10.4e9, y: 1_003e6, shear: 530e6),
            edge: .clamped
        )
        let found = `Acoustics.Plate Tests`.frequencies(try Acoustics.Plate.Modes(plate, count: 6))
        // Table 9: (0,0), (0,1), (1,0), (1,1), (0,2), (1,2).
        let measured = [57.0, 100, 130, 162, 170, 224]
        for (f, m) in zip(found, measured) {
            #expect(abs(f - m) / m < 0.03, "computed \(f) Hz, measured \(m) Hz")
        }
    }
}
