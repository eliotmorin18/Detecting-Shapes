using Test
using LinearAlgebra

# Charge la fonction depuis src/filtration.jl
# Chemin robuste depuis le dossier test/
include(joinpath(@__DIR__, "..", "src", "filtration.jl"))

@testset "pairwise_distances - cas simple" begin
    # 3 points en 3D
    points = [
        0.0 0.0 0.0;
        1.0 0.0 0.0;
        0.0 1.0 0.0
    ]

    D = pairwise_distances(points)

    # Attendu :
    # d(p1,p2)=1, d(p1,p3)=1, d(p2,p3)=sqrt(2)
    expected = [
        0.0 1.0 1.0;
        1.0 0.0 sqrt(2);
        1.0 sqrt(2) 0.0
    ]

    @test size(D) == (3, 3)
    @test D ≈ expected
    @test D ≈ D'                 # symétrie
    @test all(diag(D) .== 0.0)   # diagonale nulle
end
