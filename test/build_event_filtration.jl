using Test
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "filtration.jl"))

@testset "build_event_filtration - 3 points" begin
    points = [
        0.0 0.0 0.0;
        1.0 0.0 0.0;
        0.0 1.0 0.0
    ]
    D = pairwise_distances(points)

    @testset "eps_max = 2.0 (tout inclus)" begin
        events = build_event_filtration(D; eps_max=2.0)

        # 3 sommets + 3 arêtes + 1 triangle = 7
        @test length(events) == 7

        # Vérifier présence des sommets
        verts0 = [ev.verts for ev in events if ev.dim == 0]
        @test sort(verts0) == [(1,0,0), (2,0,0), (3,0,0)]
        @test all(ev.eps == 0.0 for ev in events if ev.dim == 0)

        # Vérifier arêtes (eps = distance)
        edges = Dict(ev.verts => ev.eps for ev in events if ev.dim == 1)
        @test edges[(1,2,0)] ≈ 1.0
        @test edges[(1,3,0)] ≈ 1.0
        @test edges[(2,3,0)] ≈ sqrt(2)

        # Vérifier triangle (eps = max des 3 arêtes)
        tris = [ev for ev in events if ev.dim == 2]
        @test length(tris) == 1
        @test tris[1].verts == (1,2,3)
        @test tris[1].eps ≈ sqrt(2)

        # Vérifier tri global: eps croissant puis dim croissant
        keys = [(ev.eps, ev.dim, ev.verts...) for ev in events]
        @test keys == sort(keys)

        # Vérifier que tout est borné par eps_max
        @test all(ev.eps <= 2.0 + 1e-12 for ev in events)
    end

    @testset "eps_max = 1.0 (triangle exclu, arête (2,3) exclue)" begin
        events = build_event_filtration(D; eps_max=1.0)

        # 3 sommets + 2 arêtes (12 et 13) = 5
        @test length(events) == 5

        edges = sort([ev.verts for ev in events if ev.dim == 1])
        @test edges == [(1,2,0), (1,3,0)]

        @test isempty([ev for ev in events if ev.dim == 2])  # pas de triangle

        # tri toujours valide
        keys = [(ev.eps, ev.dim, ev.verts...) for ev in events]
        @test keys == sort(keys)

        @test all(ev.eps <= 1.0 + 1e-12 for ev in events)
    end
end
