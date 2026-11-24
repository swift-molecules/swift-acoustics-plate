extension Acoustics.Plate {
    /// How a plate's outline is held: the edges of the mesh that only one face uses.
    public enum Edge: Hashable, Sendable, CaseIterable {
        /// Nothing holds it.
        case free
        /// Held from moving out of the plane, free to turn: the soft simple support.
        case pinned
        /// Held from moving out of the plane and from turning about the edge's normal, free to turn about the edge:
        /// the hard simple support, whose modes on a rectangle are Navier's.
        case supported
        /// Held from moving and from turning.
        case clamped
    }
}
