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
        /// How the outline is held.
        public var edge: Edge
        /// The vertices of the outline held otherwise than `edge`, by their index in the mesh's vertices: a
        /// cantilever is a free edge with the vertices of one side clamped.
        public var exceptions: [Int: Edge]

        @inlinable
        public init(
            mesh: Mesh<Point<2, Double>>,
            thickness: Measure<1, Space.Physical>.Value<Double>,
            material: Material,
            edge: Edge,
            exceptions: [Int: Edge] = [:]
        ) {
            self.mesh = mesh
            self.thickness = thickness
            self.material = material
            self.edge = edge
            self.exceptions = exceptions
        }
    }
}

extension Acoustics.Plate: Hashable {}
extension Acoustics.Plate: Sendable {}
