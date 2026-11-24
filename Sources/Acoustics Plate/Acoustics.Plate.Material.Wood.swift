public import Density
public import Modulus
public import Poisson

extension Acoustics.Plate.Material {
    /// The elastic constants of a clear wood at 12% moisture content along its three axes: longitudinal (along the
    /// grain), radial (across the growth rings) and tangential (along them).
    public struct Wood: Hashable, Sendable {
        /// The mass of a cubic metre at 12% moisture content.
        public var density: Density.Value<Double>
        /// E_L, E_R and E_T.
        public var young: Young
        /// G_LR, G_LT and G_RT.
        public var shear: Shear
        /// μ_LR and μ_LT, the narrowing across the grain per stretch along it.
        public var poisson: Poisson

        @inlinable
        public init(density: Density.Value<Double>, young: Young, shear: Shear, poisson: Poisson) {
            self.density = density
            self.young = young
            self.shear = shear
            self.poisson = poisson
        }
    }
}

extension Acoustics.Plate.Material.Wood {
    /// How a board is cut from the log, which decides the wood's axes that lie in the plate.
    public enum Cut: Hashable, Sendable, CaseIterable {
        /// The plate lies in the longitudinal–radial plane, as a soundboard does.
        case quarter
        /// The plate lies in the longitudinal–tangential plane.
        case flat
    }

    public struct Young: Hashable, Sendable {
        public var longitudinal: Modulus.Young.Value<Double>
        public var radial: Modulus.Young.Value<Double>
        public var tangential: Modulus.Young.Value<Double>

        @inlinable
        public init(
            longitudinal: Modulus.Young.Value<Double>,
            radial: Modulus.Young.Value<Double>,
            tangential: Modulus.Young.Value<Double>
        ) {
            self.longitudinal = longitudinal
            self.radial = radial
            self.tangential = tangential
        }
    }

    public struct Shear: Hashable, Sendable {
        public var longitudinalRadial: Modulus.Shear.Value<Double>
        public var longitudinalTangential: Modulus.Shear.Value<Double>
        public var radialTangential: Modulus.Shear.Value<Double>

        @inlinable
        public init(
            longitudinalRadial: Modulus.Shear.Value<Double>,
            longitudinalTangential: Modulus.Shear.Value<Double>,
            radialTangential: Modulus.Shear.Value<Double>
        ) {
            self.longitudinalRadial = longitudinalRadial
            self.longitudinalTangential = longitudinalTangential
            self.radialTangential = radialTangential
        }
    }

    public struct Poisson: Hashable, Sendable {
        public var longitudinalRadial: Poisson::Poisson.Value<Double>
        public var longitudinalTangential: Poisson::Poisson.Value<Double>

        @inlinable
        public init(
            longitudinalRadial: Poisson::Poisson.Value<Double>, longitudinalTangential: Poisson::Poisson.Value<Double>
        ) {
            self.longitudinalRadial = longitudinalRadial
            self.longitudinalTangential = longitudinalTangential
        }
    }
}

extension Acoustics.Plate.Material.Wood {
    /// The wood the USDA Wood Handbook's tables give, at 12% moisture content.
    ///
    /// - `specificGravity`: Table 5–3a, on ovendry mass and volume at 12% moisture content.
    /// - `modulus`: Table 5–3a's modulus of elasticity in megapascals, from static bending. E_L is it increased by
    ///   10%, as footnote a of Table 5–1 and footnote c of Table 5–3a give, to remove the effect of shear deflection.
    /// - `ratios`: Table 5–1's E_T/E_L, E_R/E_L, G_LR/E_L, G_LT/E_L and G_RT/E_L.
    /// - `poisson`: Table 5–2's μ_LR and μ_LT.
    ///
    /// The density is ρ_w G₁₂ (1 + 12/100) with ρ_w = 1,000 kg m⁻³, Equation 4–14.
    @inlinable
    init(
        specificGravity: Double,
        modulus: Double,
        ratios: (tangential: Double, radial: Double, lr: Double, lt: Double, rt: Double),
        poisson: (lr: Double, lt: Double)
    ) {
        let longitudinal = modulus * 1_000_000 * 1.1
        self.init(
            density: Density.amount(1_000 * specificGravity * 1.12),
            young: Young(
                longitudinal: Modulus.Young.pascals(longitudinal),
                radial: Modulus.Young.pascals(longitudinal * ratios.radial),
                tangential: Modulus.Young.pascals(longitudinal * ratios.tangential)
            ),
            shear: Shear(
                longitudinalRadial: Modulus.Shear.pascals(longitudinal * ratios.lr),
                longitudinalTangential: Modulus.Shear.pascals(longitudinal * ratios.lt),
                radialTangential: Modulus.Shear.pascals(longitudinal * ratios.rt)
            ),
            poisson: Poisson(
                longitudinalRadial: Poisson::Poisson.amount(poisson.lr),
                longitudinalTangential: Poisson::Poisson.amount(poisson.lt)
            )
        )
    }
}

// Every value below is from Forest Products Laboratory (2021), Wood handbook—wood as an engineering material,
// General Technical Report FPL-GTR-282, Madison, WI: USDA Forest Service, Forest Products Laboratory, chapter 5:
// Table 5–1 (page 5–2), Table 5–2 (page 5–3) and the 12% rows of Table 5–3a (pages 5–6 to 5–8).
// Only species whose rows in Tables 5–1 and 5–2 are complete for a plate are given.

extension Acoustics.Plate.Material.Wood {
    /// Sitka spruce (Picea sitchensis), the soundboard wood. Table 5–3a page 5–8; Table 5–1 "Spruce, Sitka";
    /// Table 5–2 "Spruce, Sitka".
    public static var sitkaSpruce: Self {
        Self(
            specificGravity: 0.40,
            modulus: 10_800,
            ratios: (tangential: 0.043, radial: 0.078, lr: 0.064, lt: 0.061, rt: 0.003),
            poisson: (lr: 0.372, lt: 0.467)
        )
    }

    /// Engelmann spruce (Picea engelmannii). Table 5–3a page 5–8; Tables 5–1 and 5–2 "Spruce, Engelmann".
    public static var engelmannSpruce: Self {
        Self(
            specificGravity: 0.35,
            modulus: 8_900,
            ratios: (tangential: 0.059, radial: 0.128, lr: 0.124, lt: 0.120, rt: 0.010),
            poisson: (lr: 0.422, lt: 0.462)
        )
    }

    /// Western redcedar (Thuja plicata). Table 5–3a page 5–7 "Western redcedar"; Tables 5–1 and 5–2
    /// "Cedar, western red".
    public static var westernRedcedar: Self {
        Self(
            specificGravity: 0.32,
            modulus: 7_700,
            ratios: (tangential: 0.055, radial: 0.081, lr: 0.087, lt: 0.086, rt: 0.005),
            poisson: (lr: 0.378, lt: 0.296)
        )
    }

    /// Douglas-fir (Pseudotsuga menziesii), Coast. Table 5–3a page 5–7 "Douglas-fir, Coast"; Tables 5–1 and 5–2
    /// "Douglas-fir".
    public static var douglasFir: Self {
        Self(
            specificGravity: 0.48,
            modulus: 13_400,
            ratios: (tangential: 0.050, radial: 0.068, lr: 0.064, lt: 0.078, rt: 0.007),
            poisson: (lr: 0.292, lt: 0.449)
        )
    }

    /// Yellow-poplar (Liriodendron tulipifera). Table 5–3a page 5–6; Tables 5–1 and 5–2 "Yellow-poplar".
    public static var yellowPoplar: Self {
        Self(
            specificGravity: 0.42,
            modulus: 10_900,
            ratios: (tangential: 0.043, radial: 0.092, lr: 0.075, lt: 0.069, rt: 0.011),
            poisson: (lr: 0.318, lt: 0.392)
        )
    }

    /// Black walnut (Juglans nigra). Table 5–3a page 5–6; Tables 5–1 and 5–2 "Walnut, black".
    public static var blackWalnut: Self {
        Self(
            specificGravity: 0.55,
            modulus: 11_600,
            ratios: (tangential: 0.056, radial: 0.106, lr: 0.085, lt: 0.062, rt: 0.021),
            poisson: (lr: 0.495, lt: 0.632)
        )
    }
}
