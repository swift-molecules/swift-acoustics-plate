@testable import Acoustics_Plate
import Testing

extension `Acoustics.Plate Tests` {
    /// Plates whose modes are known in closed form.
    @Suite struct Closed {}
}

extension `Acoustics.Plate Tests`.Closed {
    /// The lowest natural frequency of the (m, n) mode of a hard simply supported orthotropic Mindlin plate,
    /// a × b: Navier's solution w = W sin αx sin βy, βx = X cos αx sin βy, βy = Y sin αx cos βy, α = mπ/a,
    /// β = nπ/b, whose amplitudes solve a 3 × 3 pencil, solved here by swift-eigen.
    static func navier(_ plate: Acoustics.Plate, a: Double, b: Double, m: Int, n: Int) throws -> Double {
        let material = plate.material
        let h = plate.thickness.underlying
        let ex = material.young.x.underlying
        let ey = material.young.y.underlying
        let nxy = material.poisson.underlying
        let denominator = 1 - nxy * nxy * ey / ex
        let bend = h * h * h / 12
        let d11 = bend * ex / denominator
        let d22 = bend * ey / denominator
        let d12 = bend * nxy * ey / denominator
        let d66 = bend * material.shear.xy.underlying
        let s1 = 5.0 / 6 * h * material.shear.xz.underlying
        let s2 = 5.0 / 6 * h * material.shear.yz.underlying
        let alpha = Double(m) * .pi / a
        let beta = Double(n) * .pi / b
        let rho = material.density.underlying
        typealias Entry = Sparse.Symmetric<Double>.Entry
        let stiffness = try Sparse.Symmetric<Double>(
            order: 3,
            entries: [
                Entry(0, 0, s1 * alpha * alpha + s2 * beta * beta),
                Entry(0, 1, s1 * alpha),
                Entry(0, 2, s2 * beta),
                Entry(1, 1, d11 * alpha * alpha + d66 * beta * beta + s1),
                Entry(1, 2, (d12 + d66) * alpha * beta),
                Entry(2, 2, d22 * beta * beta + d66 * alpha * alpha + s2),
            ]
        )
        let mass = try Sparse.Symmetric<Double>(
            order: 3,
            entries: [Entry(0, 0, rho * h), Entry(1, 1, rho * h * h * h / 12), Entry(2, 2, rho * h * h * h / 12)]
        )
        let lowest = try Eigen.Pencil(stiffness, mass).lowest(1)[0].value
        return lowest.squareRoot() / (2 * .pi)
    }

    @Test
    func `A thin square simply supported plate rings at Navier's frequencies, the repeated ones twice`() throws {
        let a = 0.5
        let h = 0.001
        let plate = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.rectangle(a, a, columns: 32, rows: 32),
            thickness: Space.Physical.metres(h),
            material: `Acoustics.Plate Tests`.aluminium,
            edge: .supported
        )
        let modes = try Acoustics.Plate.Modes(plate, count: 6)
        let found = `Acoustics.Plate Tests`.frequencies(modes)
        // Kirchhoff: f = (π / 2) (m² + n²) / a² √(D / ρh), for (1,1), (1,2), (2,1), (2,2), (1,3), (3,1).
        let c = `Acoustics.Plate Tests`.flexural(plate.material, thickness: h) * .pi / 2 / (a * a)
        let expected = [2.0, 5, 5, 8, 10, 10].map { c * $0 }
        for (f, e) in zip(found, expected) {
            #expect(abs(f - e) / e < 0.01, "found \(f) Hz, Navier \(e) Hz")
        }
    }

    @Test
    func `A quartersawn Sitka spruce panel rings at the Mindlin–Navier frequencies of its first ten modes`() throws {
        let a = 0.40
        let b = 0.25
        let plate = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.rectangle(a, b, columns: 80, rows: 50),
            thickness: Space.Physical.metres(0.003),
            material: Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter),
            edge: .supported
        )
        var expected: [Double] = []
        for m in 1...8 {
            for n in 1...6 { expected.append(try Self.navier(plate, a: a, b: b, m: m, n: n)) }
        }
        expected.sort()
        let found = `Acoustics.Plate Tests`.frequencies(try Acoustics.Plate.Modes(plate, count: 10))
        for (f, e) in zip(found, expected) {
            #expect(abs(f - e) / e < 0.01, "found \(f) Hz, Navier \(e) Hz")
        }
    }

    @Test
    func `A clamped circular plate rings at Leissa's frequencies`() throws {
        let radius = 0.1
        let h = 0.001
        let plate = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.disk(radius: radius, cells: 16),
            thickness: Space.Physical.metres(h),
            material: `Acoustics.Plate Tests`.aluminium,
            edge: .clamped
        )
        let found = `Acoustics.Plate Tests`.frequencies(try Acoustics.Plate.Modes(plate, count: 6))
        // Leissa (1969), Vibration of Plates, NASA SP-160, Table 2.1, page 8: ω a² √(ρh / D) for (n, s) = (0, 0), (1, 0)
        // twice, (2, 0) twice, (0, 1).
        let lambda = [10.2158, 21.26, 21.26, 34.88, 34.88, 39.771]
        let c = `Acoustics.Plate Tests`.flexural(plate.material, thickness: h) / (radius * radius) / (2 * .pi)
        for (f, l) in zip(found, lambda) {
            #expect(abs(f - c * l) / (c * l) < 0.01, "found \(f) Hz, Leissa \(c * l) Hz")
        }
    }

    @Test
    func `A free square plate has three rigid motions and then Leissa's frequencies`() throws {
        let a = 0.3
        let h = 0.001
        let plate = Acoustics.Plate(
            mesh: `Acoustics.Plate Tests`.rectangle(a, a, columns: 24, rows: 24),
            thickness: Space.Physical.metres(h),
            material: `Acoustics.Plate Tests`.aluminium,
            edge: .free
        )
        let found = `Acoustics.Plate Tests`.frequencies(try Acoustics.Plate.Modes(plate, count: 6))
        let c = `Acoustics.Plate Tests`.flexural(plate.material, thickness: h) / (a * a) / (2 * .pi)
        for f in found.prefix(3) { #expect(f < 1e-3 * c) }
        // Leissa (1969), Vibration of Plates, NASA SP-160, Table 4.65, pages 104 to 106, the completely free
        // square with ν = 0.3: ω a² √(ρh / D) = 13.4728, 19.5961, 24.2702.
        let lambda = [13.4728, 19.5961, 24.2702]
        for (f, l) in zip(found.dropFirst(3), lambda) {
            #expect(abs(f - c * l) / (c * l) < 0.005, "found \(f) Hz, Leissa \(c * l) Hz")
        }
    }

    @Test
    func `The fundamental's error falls by more than three each time the mesh is halved`() throws {
        let a = 0.5
        let h = 0.001
        let exact = `Acoustics.Plate Tests`.flexural(`Acoustics.Plate Tests`.aluminium, thickness: h) * .pi / (a * a)
        var errors: [Double] = []
        for cells in [4, 8, 16, 32] {
            let plate = Acoustics.Plate(
                mesh: `Acoustics.Plate Tests`.rectangle(a, a, columns: cells, rows: cells),
                thickness: Space.Physical.metres(h),
                material: `Acoustics.Plate Tests`.aluminium,
                edge: .supported
            )
            let f = try Acoustics.Plate.Modes(plate, count: 1)[0].frequency.underlying
            errors.append((f - exact) / exact)
        }
        for (coarse, fine) in zip(errors, errors.dropFirst()) {
            #expect(abs(fine) < abs(coarse) / 3, "relative errors \(errors)")
        }
        #expect(abs(errors.last!) < 0.002)
    }
}
