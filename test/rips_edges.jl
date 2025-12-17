using Test
using LinearAlgebra

# Charge les fonctions UNE SEULE FOIS (en haut)
include(joinpath(@__DIR__, "..", "src", "filtration.jl"))

@testset "rips_edges" begin
    @testset "cas simple" begin
        points = [
            0.0 0.0 0.0;
            1.0 0.0 0.0;
            0.0 1.0 0.0
        ]

        D = pairwise_distances(points)

        edges1 = rips_edges(D, 1.0)
        @test sort(edges1) == [(1,2), (1,3)]

        edges2 = rips_edges(D, 1.5)
        @test sort(edges2) == [(1,2), (1,3), (2,3)]

        edges3 = rips_edges(D, 0.9)
        @test isempty(edges3)
    end

    @testset "cas N=2" begin
        D2 = [0.0 0.5; 0.5 0.0]
        @test rips_edges(D2, 0.4) == Tuple{Int,Int}[]
        @test rips_edges(D2, 0.5) == [(1,2)]
    end
end
