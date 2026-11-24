// swift-tools-version: 6.4

import PackageDescription

let package = Package(
    name: "swift-acoustics-plate",
    platforms: [
        .macOS(.v27),
        .iOS(.v27),
        .tvOS(.v27),
        .watchOS(.v27),
        .visionOS(.v27),
    ],
    products: [
        .library(name: "Acoustics Plate", targets: ["Acoustics Plate"])
    ],
    dependencies: [
        .package(url: "https://github.com/swift-molecules/swift-acoustics.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-sparse.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-eigen.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-mesh.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-numeric.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-density.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-modulus.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-poisson.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-frequency.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-bessel.git", branch: "main"),
    ],
    targets: [
        .target(
            name: "Acoustics Plate",
            dependencies: [
                .product(name: "Acoustics", package: "swift-acoustics"),
                .product(name: "Sparse", package: "swift-sparse"),
                .product(name: "Eigen", package: "swift-eigen"),
                .product(name: "Mesh", package: "swift-mesh"),
                .product(name: "Numeric", package: "swift-numeric"),
                .product(name: "Density", package: "swift-density"),
                .product(name: "Modulus", package: "swift-modulus"),
                .product(name: "Poisson", package: "swift-poisson"),
                .product(name: "Frequency", package: "swift-frequency"),
            ]
        ),
        .testTarget(
            name: "Acoustics Plate Tests",
            dependencies: [
                .target(name: "Acoustics Plate"),
                .product(name: "Bessel", package: "swift-bessel"),
                .product(name: "Numeric", package: "swift-numeric"),
            ]
        ),
    ],
    swiftLanguageModes: [.v6]
)

for target in package.targets where ![.system, .binary, .plugin, .macro].contains(target.type) {
    target.swiftSettings =
        (target.swiftSettings ?? []) + [
            .strictMemorySafety(),
            .enableUpcomingFeature("ExistentialAny"),
            .enableUpcomingFeature("InternalImportsByDefault"),
            .enableUpcomingFeature("MemberImportVisibility"),
            .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
            .enableUpcomingFeature("InferIsolatedConformances"),
            .enableExperimentalFeature("Lifetimes"),
            .treatAllWarnings(as: .error),
        ]
}
