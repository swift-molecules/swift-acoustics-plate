extension Acoustics.Plate {
    /// One quadrilateral of Reissner–Mindlin plate: deflection w and the rotations βx, βy of the normal at each
    /// corner, so that a point at height z moves z βx along x and z βy along y. Bending is integrated at 2 × 2 Gauss
    /// points; the transverse shear is MITC4's (Bathe and Dvorkin 1985), its covariant strains tied at the edges'
    /// midpoints so that the element does not lock when the plate is thin. The mass is consistent and includes the
    /// normal's rotary inertia ρh³/12.
    @usableFromInline
    internal struct Element {
        /// Corner coordinates, counterclockwise.
        @usableFromInline let x: [Double]
        @usableFromInline let y: [Double]

        @usableFromInline
        init(x: [Double], y: [Double]) {
            self.x = x
            self.y = y
        }
    }
}

extension Acoustics.Plate.Element {
    /// The corners' signs in the reference square.
    @usableFromInline static let corners: [(Double, Double)] = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
    /// The 2 × 2 Gauss points, each of weight one.
    @usableFromInline static let gauss: [(Double, Double)] = {
        let g = 1 / Double(3).squareRoot()
        return [(-g, -g), (g, -g), (g, g), (-g, g)]
    }()

    /// The shape functions and their derivatives in the reference square at (ξ, η).
    @usableFromInline
    static func shape(_ xi: Double, _ eta: Double) -> (n: [Double], dxi: [Double], deta: [Double]) {
        var n = [Double](repeating: 0, count: 4)
        var dxi = n
        var deta = n
        for i in 0..<4 {
            let (a, b) = corners[i]
            n[i] = 0.25 * (1 + a * xi) * (1 + b * eta)
            dxi[i] = 0.25 * a * (1 + b * eta)
            deta[i] = 0.25 * b * (1 + a * xi)
        }
        return (n, dxi, deta)
    }

    /// The Jacobian [[x,ξ, y,ξ], [x,η, y,η]] at (ξ, η).
    @usableFromInline
    func jacobian(_ dxi: [Double], _ deta: [Double]) -> (xXi: Double, yXi: Double, xEta: Double, yEta: Double) {
        var j = (xXi: 0.0, yXi: 0.0, xEta: 0.0, yEta: 0.0)
        for i in 0..<4 {
            j.xXi += dxi[i] * x[i]
            j.yXi += dxi[i] * y[i]
            j.xEta += deta[i] * x[i]
            j.yEta += deta[i] * y[i]
        }
        return j
    }

    /// The smallest Jacobian determinant over the Gauss points and the corners: positive for a counterclockwise
    /// convex quadrilateral.
    @usableFromInline
    var determinant: Double {
        var least = Double.infinity
        for (xi, eta) in Self.gauss + Self.corners {
            let s = Self.shape(xi, eta)
            let j = jacobian(s.dxi, s.deta)
            least = Swift.min(least, j.xXi * j.yEta - j.yXi * j.xEta)
        }
        return least
    }

    /// The covariant shear strain γ_ξz (`along` true) or γ_ηz at (ξ, η) as a row over the twelve unknowns.
    @usableFromInline
    func covariant(_ xi: Double, _ eta: Double, along: Bool) -> [Double] {
        let s = Self.shape(xi, eta)
        let j = jacobian(s.dxi, s.deta)
        let d = along ? s.dxi : s.deta
        let dx = along ? j.xXi : j.xEta
        let dy = along ? j.yXi : j.yEta
        var row = [Double](repeating: 0, count: 12)
        for i in 0..<4 {
            row[3 * i] = d[i]
            row[3 * i + 1] = s.n[i] * dx
            row[3 * i + 2] = s.n[i] * dy
        }
        return row
    }

    /// The 12 × 12 stiffness and mass matrices, row-major, for a plate of thickness `h`, density `rho`, reduced
    /// in-plane stiffness `q` = (Q11, Q12, Q22, Q66) and transverse shear stiffness `s` = (κ G_xz, κ G_yz).
    @usableFromInline
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
        for (xi, eta) in Self.gauss {
            let sh = Self.shape(xi, eta)
            let j = jacobian(sh.dxi, sh.deta)
            let det = j.xXi * j.yEta - j.yXi * j.xEta
            // J⁻¹ = [[y,η, −y,ξ], [−x,η, x,ξ]] / det maps (∂/∂ξ, ∂/∂η) to (∂/∂x, ∂/∂y).
            let inverse = (a: j.yEta / det, b: -j.yXi / det, c: -j.xEta / det, d: j.xXi / det)
            var dx = [Double](repeating: 0, count: 4)
            var dy = dx
            for i in 0..<4 {
                dx[i] = inverse.a * sh.dxi[i] + inverse.b * sh.deta[i]
                dy[i] = inverse.c * sh.dxi[i] + inverse.d * sh.deta[i]
            }
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
                gxz[a] = inverse.a * gxi[a] + inverse.b * geta[a]
                gyz[a] = inverse.c * gxi[a] + inverse.d * geta[a]
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
                    let nn = det * sh.n[i] * sh.n[l]
                    m[12 * (3 * i) + 3 * l] += translational * nn
                    m[12 * (3 * i + 1) + 3 * l + 1] += rotary * nn
                    m[12 * (3 * i + 2) + 3 * l + 2] += rotary * nn
                }
            }
        }
        return (k, m)
    }

    /// ∫ Nᵢ dA for each corner, ∫ x dA and ∫ y dA, and the area.
    @usableFromInline
    var moments: (weights: [Double], x: Double, y: Double, area: Double) {
        var weights = [Double](repeating: 0, count: 4)
        var mx = 0.0
        var my = 0.0
        var area = 0.0
        for (xi, eta) in Self.gauss {
            let sh = Self.shape(xi, eta)
            let j = jacobian(sh.dxi, sh.deta)
            let det = j.xXi * j.yEta - j.yXi * j.xEta
            for i in 0..<4 {
                weights[i] += det * sh.n[i]
                mx += det * sh.n[i] * x[i]
                my += det * sh.n[i] * y[i]
            }
            area += det
        }
        return (weights, mx, my, area)
    }

    /// ∫ Nᵢ Nⱼ dA, row-major 4 × 4.
    @usableFromInline
    var overlap: [Double] {
        var result = [Double](repeating: 0, count: 16)
        for (xi, eta) in Self.gauss {
            let sh = Self.shape(xi, eta)
            let j = jacobian(sh.dxi, sh.deta)
            let det = j.xXi * j.yEta - j.yXi * j.xEta
            for i in 0..<4 {
                for l in 0..<4 { result[4 * i + l] += det * sh.n[i] * sh.n[l] }
            }
        }
        return result
    }
}
