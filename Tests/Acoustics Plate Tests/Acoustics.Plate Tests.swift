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

    /// A disk of radius `radius` about the origin: a square of half-side 0.4 radius in the middle, cut into
    /// `cells` × `cells` quads, and four patches of `cells` × `cells` quads between its sides and the circle, each
    /// point of a side projected outward onto the circle.
    static func disk(radius: Double, cells: Int) -> Mesh<Point<2, Double>> {
        let s = 0.4 * radius
        let count = Cardinal(UInt(cells))
        func parameter(_ index: Ordinal) -> Double { Double(2 * Int(bitPattern: index) - cells) / Double(cells) }
        var mesh = Mesh(columns: count, rows: count) { column, row in
            Point<2, Double>(x: s * parameter(column), y: s * parameter(row))
        }
        let sides: [(Double) -> (Double, Double)] = [
            { t in (s, s * t) },
            { t in (-s * t, s) },
            { t in (-s, -s * t) },
            { t in (s * t, -s) },
        ]
        for side in sides {
            mesh = mesh.appending(
                Mesh(columns: count, rows: count) { column, row in
                    let (x, y) = side(parameter(column))
                    let r = Double(Int(bitPattern: row)) / Double(cells)
                    let scale = radius / (x * x + y * y).squareRoot()
                    return Point<2, Double>(x: (1 - r) * x + r * scale * x, y: (1 - r) * y + r * scale * y)
                }
            )
        }
        return mesh.welded
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
