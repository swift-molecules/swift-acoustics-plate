public import Density
public import Modulus
public import Poisson

extension Acoustics.Plate {
    /// An orthotropic material whose axes of symmetry are the plate's: x and y in its plane, z through its thickness.
    ///
    /// An isotropic material is the case where both Young's moduli are equal and every shear modulus is
    /// E / 2(1 + ν).
    public struct Material {
        /// The mass of a cubic metre.
        public var density: Density.Value<Double>
        /// Young's moduli along x and y.
        public var young: Young
        /// The shear moduli in the x–y, x–z and y–z planes.
        public var shear: Shear
        /// Poisson's ratio ν_xy: the narrowing along y per stretch along x. ν_yx follows as ν_xy E_y / E_x.
        public var poisson: Poisson.Value<Double>

        @inlinable
        public init(
            density: Density.Value<Double>,
            young: Young,
            shear: Shear,
            poisson: Poisson.Value<Double>
        ) {
            self.density = density
            self.young = young
            self.shear = shear
            self.poisson = poisson
        }
    }
}

extension Acoustics.Plate.Material {
    /// Young's moduli along the plate's two axes.
    public struct Young: Hashable, Sendable {
        public var x: Modulus.Young.Value<Double>
        public var y: Modulus.Young.Value<Double>

        @inlinable
        public init(x: Modulus.Young.Value<Double>, y: Modulus.Young.Value<Double>) {
            self.x = x
            self.y = y
        }
    }

    /// The shear moduli in the plate's three planes.
    public struct Shear: Hashable, Sendable {
        /// In the plate's plane.
        public var xy: Modulus.Shear.Value<Double>
        /// Through the thickness, along x.
        public var xz: Modulus.Shear.Value<Double>
        /// Through the thickness, along y.
        public var yz: Modulus.Shear.Value<Double>

        @inlinable
        public init(
            xy: Modulus.Shear.Value<Double>,
            xz: Modulus.Shear.Value<Double>,
            yz: Modulus.Shear.Value<Double>
        ) {
            self.xy = xy
            self.xz = xz
            self.yz = yz
        }
    }
}

extension Acoustics.Plate.Material: Hashable {}
extension Acoustics.Plate.Material: Sendable {}

extension Acoustics.Plate.Material {
    /// The isotropic material of this density, Young's modulus and Poisson's ratio, every shear modulus E / 2(1 + ν).
    @inlinable
    public static func isotropic(
        density: Density.Value<Double>,
        young: Modulus.Young.Value<Double>,
        poisson: Poisson.Value<Double>
    ) -> Self {
        let shear = Modulus.Shear.pascals(young.underlying / (2 * (1 + poisson.underlying)))
        return Self(
            density: density,
            young: Young(x: young, y: young),
            shear: Shear(xy: shear, xz: shear, yz: shear),
            poisson: poisson
        )
    }

    /// This wood, cut so that its grain runs along x.
    ///
    /// Quartersawn, the plate lies in the wood's longitudinal–radial plane and its thickness is tangential:
    /// E_x = E_L, E_y = E_R, G_xy = G_LR, G_xz = G_LT, G_yz = G_RT, ν_xy = μ_LR. Flatsawn, it lies in the
    /// longitudinal–tangential plane and its thickness is radial: E_y = E_T, G_xy = G_LT, G_xz = G_LR,
    /// G_yz = G_RT, ν_xy = μ_LT.
    @inlinable
    public init(_ wood: Wood, cut: Wood.Cut) {
        switch cut {
        case .quarter:
            self.init(
                density: wood.density,
                young: Young(x: wood.young.longitudinal, y: wood.young.radial),
                shear: Shear(xy: wood.shear.longitudinalRadial, xz: wood.shear.longitudinalTangential, yz: wood.shear.radialTangential),
                poisson: wood.poisson.longitudinalRadial
            )
        case .flat:
            self.init(
                density: wood.density,
                young: Young(x: wood.young.longitudinal, y: wood.young.tangential),
                shear: Shear(xy: wood.shear.longitudinalTangential, xz: wood.shear.longitudinalRadial, yz: wood.shear.radialTangential),
                poisson: wood.poisson.longitudinalTangential
            )
        }
    }
}
