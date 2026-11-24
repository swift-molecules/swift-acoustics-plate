public import Density
public import Frequency

extension Acoustics.Plate {
    /// One mode of a plate: how fast it rings, its shape, the mass that moves with it, and how well it radiates.
    public struct Mode {
        /// The natural frequency.
        public let frequency: Frequency.Value<Double>
        /// The deflection at each vertex of the plate's mesh, in its order, scaled so that the largest is one and
        /// positive. A held vertex has zero.
        public let shape: [Double]
        /// The modal mass at that scale: the kinetic energy of the mode is one half this mass times the square of
        /// the velocity of its largest deflection. A resonator mode with this mass, struck where the shape is s,
        /// moves the plate there s times as far as its own displacement.
        public let mass: Mass.Value<Double>
        /// The signed area the shape sweeps, ∫ shape dA: the volume the plate pushes out per metre of its largest
        /// deflection, which sets how strongly the mode radiates at low frequency.
        public let area: Area<Space.Physical>.Value<Double>
        /// The radiation efficiency at the natural frequency, set in an infinite rigid baffle into air on one side:
        /// the power radiated over ρ₀ c times the plate's area times the mean square normal velocity.
        public let efficiency: Double

        @inlinable
        public init(
            frequency: Frequency.Value<Double>,
            shape: [Double],
            mass: Mass.Value<Double>,
            area: Area<Space.Physical>.Value<Double>,
            efficiency: Double
        ) {
            self.frequency = frequency
            self.shape = shape
            self.mass = mass
            self.area = area
            self.efficiency = efficiency
        }
    }
}

extension Acoustics.Plate.Mode: Hashable {}
extension Acoustics.Plate.Mode: Sendable {}
