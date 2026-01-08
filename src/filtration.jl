# src/filtration.jl
module Filtration

export rips_filtration, cech_filtration, quantile_distances

using LinearAlgebra 
using Combinatorics 
using Statistics

function rips_complex(S, r) 
    # builds the Rips complex for the point cloud S with radius r 
    
    # INPUT
    # S = point cloud 
    # r = radius 

    # EXAMPLE INPUT 
    # S = [(0, 0), (3, 9), (5, -1), (9, 4), (7, -5)]
    # r = 1 
    
    # OUTPUT

    
    # RC = Rips complex 
    n = length(S)

    # initialize and add all vertices to the Rieps Complex 
    RC = Vector{Tuple{Int64, Vararg{Int64}}}() 
    append!(RC,[(i,) for i in 1:n])
    # later L = simplices of lower dimension
    L = Vector{Tuple{Int64, Vararg{Int64}}}() 
    # Neighbours of the simplices, here of the vertices
    N = Dict(i => [] for i in RC)
    # add all those edges, which vertices are at most r apart from each other
    for i in 1:n 
        for j in (i+1):n 
            # check the distance of the vertices 
            if norm(S[i] .- S[j]) <= r 
                # add the i,j edge to the complex 
                push!(L, (i,j))
                # add the neighbours
                push!(N[(i,)],j)
                push!(N[(j,)],i)
            end
        end
    end
    append!(RC, L)
    # add higher dimensional simplices 
    dim = 2
    # new neighbours 
    N_new = Dict(face => [] for face in L)
    
    # only add simplices with dim <= dimension of the space +1 
    max_dim = length(S[1])
    while ( !isempty(L) ) & ( dim < max_dim ) 
        # increase the dimension by one
        dim += 1
        # a new simplex of dimension dim+1 should be added 
        # if all its faces with dimension dim are contained in RC
        # <=> there is one vertec s.t. each "subface" of one face build a face with this vertex 
        new_simplices = []

        for face in L 
            subfaces = [face[1:end .!= i] for i in 1:length(face)]
            intersection = intersect((N[sf] for sf in subfaces)...)

            for vx in intersection  
                new_sx = sort((face..., vx))
                if !(new_sx in new_simplices)
                    push!(new_simplices,new_sx)
                end 
                push!(N_new[face], vx)  
            end
        end 
        # update L 
        L = new_simplices
        # update N, N_new
        N = N_new
        N_new = Dict(face => [] for face in L)
        # update RC 
        append!(RC, new_simplices)
    end

    return RC
end 

function rips_filtration(S, R)
    # builds the Rips filtration for the point cloud S with radi given by the list R 
    
    # INPUT
    # S = point cloud 
    # R = list with different radii

    # EXAMPLE INPUT 
    # S = [(0, 0), (3, 9), (5, -1), (9, 4), (7, -5)]
    # R = [1,2,3,4,5] 
    
    # OUTPUT
    # RF = Ripsfiltration = Dict(R[1] => rips_complex[S,R[1]], ... )

    RF = Dict(r => rips_complex(S,r) for r in R)
    return RF 
end 

function cech_complex(S, r) 
    # builds the Chech complex for the point cloud S with radius r 
    
    # INPUT
    # S = point cloud 
    # r = radius 

    # EXAMPLE INPUT 
    # S = [(0, 0), (3, 9), (5, -1), (9, 4), (7, -5)]
    # r = 1 
    
    # OUTPUT
    # CC = Cech complex 

    n = length(S)

    # initialize and add all vertices to the Rieps Complex 
    CC = Vector{Tuple{Int64, Vararg{Int64}}}() 
    append!(CC,[(i,) for i in 1:n])
    # later L = simplices of lower dimension 
    L = Vector{Tuple{Int64, Vararg{Int64}}}() 
    # Neighbours of the simplices, here of the vertices
    N = Dict(i => [] for i in CC)
    for i in 1:n 
        for j in (i+1):n 
            # check if the balls of radius r would intersect, i.e. dist <= 2r
            if norm(S[i] .- S[j]) <= 2r 
                # add the i,j edge to the complex 
                push!(L, (i,j))
                # add the neighbours
                push!(N[(i,)],j)
                push!(N[(j,)],i)
            end
        end
    end
    append!(CC, L)

    # add higher dimensional simplices 
    dim = 2
    # only add simplices with dim <= dimension of the space +1 
    max_dim = length(S[1])
    # new neighbours 
    N_new = Dict(face => [] for face in L)
    S_collected = [collect(p) for p in S]

    while ( !isempty(L) ) & ( dim < max_dim )
        # increase the dimension by one
        dim += 1
        # a new simplex of dimension dim+1 should be added 
        # if all its faces with dimension dim are contained in RC
        # <=> there is one vertec s.t. each "subface" of one face build a face with this vertex 
        new_simplices = []

        for face in L 
            subfaces = [face[1:end .!= i] for i in 1:length(face)]
            intersection = intersect((N[sf] for sf in subfaces)...)

            for vx in intersection  
                new_sx = sort((face..., vx))
                if !(new_sx in new_simplices)
                    # compute the radius of the minimal enclosing ball
                    r_meb = minimal_enclosing_ball(S_collected[collect(new_sx)])[1]
                    if r_meb <= r 
                        push!(new_simplices,new_sx)
                    end
                else 
                    push!(N_new[face], vx)
                end 
            end
        end 
        # update L 
        L = new_simplices
        # update N, N_new
        N = N_new
        N_new = Dict(face => [] for face in L)
        # update CC 
        append!(CC, new_simplices)
    end 

    return CC
end 

function cech_filtration(S, R)
    # builds the Cech filtration for the point cloud S with radi given by the list R 
    
    # INPUT
    # S = point cloud 
    # R = list with different radi 

    # EXAMPLE INPUT 
    # S = [(0, 0), (3, 9), (5, -1), (9, 4), (7, -5)]
    # R = [1,2,3,4,5] 
    
    # OUTPUT
    # CF = Ripsfiltration = Dict(R[1] => cech_complex[S,R[1]], ... )

    CF = Dict(r => cech_complex(S,r) for r in R)
    return CF 
end 

function minimal_enclosing_ball(S)
    # computes the origin and radius of the minimal enclosing ball containing all points of S 

    # INPUT 
    # S = point cloud 

    # EXAMPLE INPUT 
    # S = [(0, 0), (3, 9), (5, -1), (9, 4), (7, -5)]
    
    # OUTPUT r,origin
    # origin = point 
    # r = radius 

    # number of S 
    n = length(S)
    if n == 2 
        # if we have two S the radius of the ball with smalles radius is the distance of the two S
        r = norm(S[1] .- S[2])/2
        origin = S[1] .+ 1/2 .* (S[2] .- S[1])
        return r,origin
    else 
        for i in 1:n
            idx = collect(1:n)
            popat!(idx,i)
            r,origin = minimal_enclosing_ball(S[idx])

            if norm(S[i] .- origin) <= r 
                return r,origin
            end 
        end 
        # otherwise the minimal_enclosing_ball is the ball through all vertices 
        
        # dimension 
        d = length(S[1])

        A = zeros(n-1, d)
        b = zeros(n-1)
        for i in 2:n
            A[i-1, :] = 2 .* (S[i] .- S[1])
            b[i-1] = norm(S[i])^2 - norm(S[1])^2
        end
        origin = A\b 
        return norm(origin - S[1]), origin
    end 
end 

function quantile_distances(S; n_quantiles = 5)
    # computes n_quantiles between 0.05 and 1 to get good values for the radius in the filtration 

    # INPUT 
    # S = pointcloud 
    # n_quantiles = integer indicating how many values to produce 

    # EXAMPLE INPUT 
    # S = [(0, 0), (3, 9), (5, -1), (9, 4), (7, -5)]

    # OUTPUT 
    # qtls = [r_1, ... r_{n_quantiles}]

    distances = [norm(collect(S[i]) - collect(S[j])) for j in 1:length(S) for i in (j+1):length(S)]
    distances = sort(distances)

    qtls = quantile(distances, range(0.05, 1; length=n_quantiles))

    return qtls
end 

end # end module