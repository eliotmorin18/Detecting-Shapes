using Test
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "filtration.jl"))

@testset "rips_triangles" begin
    @testset "3 points - triangle apparaît quand eps >= sqrt(2)" begin
        points = [
            0.0 0.0 0.0;  # 1
            1.0 0.0 0.0;  # 2
            0.0 1.0 0.0   # 3
        ]
        D = pairwise_distances(points)

        # eps trop petit: pas toutes les arêtes (2,3) manque
        T1 = rips_triangles(D, 1.0)
        @test isempty(T1)

        # eps assez grand: toutes les arêtes présentes => un triangle
        T2 = rips_triangles(D, 1.5)
        @test sort(T2) == [(1,2,3)]
    end

    @testset "4 points - eps=1.0 => aucun triangle (diagonales absentes)" begin
        points = [
            0.0 0.0 0.0;  # 1
            1.0 0.0 0.0;  # 2
            0.0 1.0 0.0;  # 3
            1.0 1.0 0.0   # 4
        ]
        D = pairwise_distances(points)

        # eps=1.0: on a les côtés du carré, mais pas les diagonales (sqrt(2))
        # donc aucun triple n'a ses 3 arêtes <= 1.0
        T = rips_triangles(D, 1.0)
        @test isempty(T)
    end
end
