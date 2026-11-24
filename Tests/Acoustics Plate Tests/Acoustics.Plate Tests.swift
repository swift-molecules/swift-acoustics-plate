import Testing

@testable import Acoustics_Plate

@Suite struct `Acoustics.Plate Tests` {}

extension `Acoustics.Plate Tests` {
    /// A rectangle `a` × `b` with its lower left corner at the origin, cut into `columns` × `rows` quads.
    static func rectangle(_ a: Double, _ b: Double, columns: Int, rows: Int) -> Mesh<Point<2, Double>> {
        Mesh(columns: Cardinal(UInt(columns)), rows: Cardinal(UInt(rows))) { column, row in
            Point<2, Double>(
                x: a * Double(Int(bitPattern: column)) / Double(columns),
                y: b * Double(Int(bitPattern: row)) / Double(rows)
            )
        }
    }

    static let aluminium = Acoustics.Plate.Material.isotropic(
        density: Density.amount(2_700),
        young: Modulus.Young.pascals(70e9),
        poisson: Poisson.amount(0.3)
    )

    /// √(D / ρh) for an isotropic plate of this material and thickness.
    static func flexural(_ material: Acoustics.Plate.Material, thickness h: Double) -> Double {
        let e = material.young.x.underlying
        let nu = material.poisson.underlying
        let d = e * h * h * h / (12 * (1 - nu * nu))
        return (d / (material.density.underlying * h)).squareRoot()
    }

    static func frequencies(_ modes: Acoustics.Plate.Modes) -> [Double] {
        modes.map(\.frequency.underlying)
    }
}
