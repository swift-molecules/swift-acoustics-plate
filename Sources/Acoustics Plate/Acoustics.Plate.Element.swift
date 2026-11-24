internal import Galerkin

extension Acoustics.Plate {
    /// One quadrilateral of Reissner–Mindlin plate: deflection w and the rotations βx, βy of the normal at each
    /// corner, so that a point at height z moves z βx along x and z βy along y. Bending is integrated at 2 × 2 Gauss
    /// points; the transverse shear is MITC4's (Bathe and Dvorkin 1985), its covariant strains tied at the edges'
    /// midpoints so that the element does not lock when the plate is thin. The mass is consistent and includes the
    /// normal's rotary inertia ρh³/12.
    internal struct Element {
        /// Corner coordinates, counterclockwise.
        let corners: [Point<2, Double>]

        init(_ corners: [Point<2, Double>]) {
            self.corners = corners
        }
    }
}

extension Acoustics.Plate.Element {
    typealias Square = Galerkin.Reference.Quadrilateral

    /// The 2 × 2 Gauss rule, exact to degree three.
    static let quadrature = Square.quadrature(degree: 3, Double.self)

    /// The shape functions and the Jacobian at (ξ, η).
    func map(_ point: Point<2, Double>) -> (shape: Galerkin.Shape<Double>, jacobian: Galerkin.Jacobian<Double>) {
        let shape = Square.shape(at: point)
        return (shape, Galerkin.Jacobian(shape, corners: corners))
    }

    /// The smallest Jacobian determinant over the Gauss points and the corners: positive for a counterclockwise
    /// convex quadrilateral.
    var determinant: Double {
        var least = Double.infinity
        for point in Self.quadrature.points + Square.corners(Double.self) {
            least = Swift.min(least, map(point).jacobian.determinant)
        }
        return least
    }

    /// The covariant shear strain γ_ξz (`along` true) or γ_ηz at (ξ, η) as a row over the twelve unknowns.
    func covariant(_ xi: Double, _ eta: Double, along: Bool) -> [Double] {
        let (s, j) = map(Point(x: xi, y: eta))
        let d = along ? s.xi : s.eta
        let dx = along ? j.xXi : j.xEta
        let dy = along ? j.yXi : j.yEta
        var row = [Double](repeating: 0, count: 12)
        for i in 0..<4 {
            row[3 * i] = d[i]
            row[3 * i + 1] = s.values[i] * dx
            row[3 * i + 2] = s.values[i] * dy
        }
        return row
    }

    /// The 12 × 12 stiffness and mass matrices, row-major, for a plate of thickness `h`, density `rho`, reduced
    /// in-plane stiffness `q` = (Q11, Q12, Q22, Q66) and transverse shear stiffness `s` = (κ G_xz, κ G_yz).
    func matrices(
        thickness h: Double,
        density rho: Double,
        q: (q11: Double, q12: Double, q22: Double, q66: Double),
        s: (xz: Double, yz: Double)
    ) -> (stiffness: [Double], mass: [Double]) {
        var k = [Double](repeating: 0, count: 144)
        var m = [Double](repeating: 0, count: 144)
        let bend = h * h * h / 12
        let d = (d11: bend * q.q11, d12: bend * q.q12, d22: bend * q.q22, d66: bend * q.q66)
        let tie = (
            a: covariant(0, 1, along: true),
            c: covariant(0, -1, along: true),
            dd: covariant(1, 0, along: false),
            b: covariant(-1, 0, along: false)
        )
        let translational = rho * h
        let rotary = rho * h * h * h / 12
        for (point, weight) in zip(Self.quadrature.points, Self.quadrature.weights) {
            let xi = point.x
            let eta = point.y
            let (sh, j) = map(point)
            let det = weight * j.determinant
            let (dx, dy) = j.gradients(sh)
            // Curvatures κxx = βx,x, κyy = βy,y, κxy = βx,y + βy,x.
            var kxx = [Double](repeating: 0, count: 12)
            var kyy = kxx
            var kxy = kxx
            for i in 0..<4 {
                kxx[3 * i + 1] = dx[i]
                kyy[3 * i + 2] = dy[i]
                kxy[3 * i + 1] = dy[i]
                kxy[3 * i + 2] = dx[i]
            }
            // The assumed covariant shear strains, then the Cartesian ones through J⁻¹.
            var gxi = [Double](repeating: 0, count: 12)
            var geta = gxi
            var gxz = gxi
            var gyz = gxi
            for a in 0..<12 {
                gxi[a] = 0.5 * (1 + eta) * tie.a[a] + 0.5 * (1 - eta) * tie.c[a]
                geta[a] = 0.5 * (1 + xi) * tie.dd[a] + 0.5 * (1 - xi) * tie.b[a]
                (gxz[a], gyz[a]) = j.cartesian(xi: gxi[a], eta: geta[a])
            }
            for a in 0..<12 {
                let mxx = d.d11 * kxx[a] + d.d12 * kyy[a]
                let myy = d.d12 * kxx[a] + d.d22 * kyy[a]
                let mxy = d.d66 * kxy[a]
                let qx = s.xz * gxz[a]
                let qy = s.yz * gyz[a]
                for b in 0..<12 {
                    k[12 * a + b] += det * (mxx * kxx[b] + myy * kyy[b] + mxy * kxy[b] + qx * gxz[b] + qy * gyz[b])
                }
            }
            for i in 0..<4 {
                for l in 0..<4 {
                    let nn = det * sh.values[i] * sh.values[l]
                    m[12 * (3 * i) + 3 * l] += translational * nn
                    m[12 * (3 * i + 1) + 3 * l + 1] += rotary * nn
                    m[12 * (3 * i + 2) + 3 * l + 2] += rotary * nn
                }
            }
        }
        return (k, m)
    }

    /// ∫ Nᵢ dA for each corner, ∫ x dA and ∫ y dA, and the area.
    var moments: (weights: [Double], x: Double, y: Double, area: Double) {
        var weights = [Double](repeating: 0, count: 4)
        var mx = 0.0
        var my = 0.0
        var area = 0.0
        for (point, weight) in zip(Self.quadrature.points, Self.quadrature.weights) {
            let (sh, j) = map(point)
            let det = weight * j.determinant
            for i in 0..<4 {
                weights[i] += det * sh.values[i]
                mx += det * sh.values[i] * corners[i].x
                my += det * sh.values[i] * corners[i].y
            }
            area += det
        }
        return (weights, mx, my, area)
    }

    /// ∫ Nᵢ Nⱼ dA, row-major 4 × 4.
    var overlap: [Double] {
        var result = [Double](repeating: 0, count: 16)
        for (point, weight) in zip(Self.quadrature.points, Self.quadrature.weights) {
            let (sh, j) = map(point)
            let det = weight * j.determinant
            for i in 0..<4 {
                for l in 0..<4 { result[4 * i + l] += det * sh.values[i] * sh.values[l] }
            }
        }
        return result
    }
}
