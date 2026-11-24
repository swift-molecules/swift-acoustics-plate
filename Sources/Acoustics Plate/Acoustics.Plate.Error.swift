public import Eigen

extension Acoustics.Plate {
    /// Why a plate's modes could not be found.
    public enum Error: Swift.Error, Hashable, Sendable {
        /// The face at this index is not a quadrilateral, or has no area.
        case face(Int)
        /// The eigen solver failed.
        case eigen(Eigen.Error)
        /// More modes were asked for than the plate has.
        case count
    }
}
