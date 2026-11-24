public import Acoustics
public import Density
public import Mesh

extension Acoustics {
    /// A flat plate: the shape of its midsurface, how thick it is, what it is made of and how its outline is held.
    ///
    /// The midsurface is a mesh of quadrilaterals in the x–y plane, its coordinates in metres, and the material's
    /// first axis runs along x. Its modes come from `Acoustics.Plate.Modes`, computed once and kept.
    public struct Plate {
        /// The midsurface: quadrilateral faces over points in metres.
        public var mesh: Mesh<Point<2, Double>>
        /// The thickness, the same everywhere.
        public var thickness: Measure<1, Space.Physical>.Value<Double>
        /// The material, its first axis along x.
        public var material: Material
        /// How the whole outline is held.
        public var edge: Edge

        @inlinable
        public init(
            mesh: Mesh<Point<2, Double>>,
            thickness: Measure<1, Space.Physical>.Value<Double>,
            material: Material,
            edge: Edge
        ) {
            self.mesh = mesh
            self.thickness = thickness
            self.material = material
            self.edge = edge
        }
    }
}

extension Acoustics.Plate: Hashable {}
extension Acoustics.Plate: Sendable {}
