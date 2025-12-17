""" functions to build the rips filtration """

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

    #initialize and add all vertices to the Rieps Complex 
    RC = Vector{Tuple{Int64, Vararg{Int64}}}() 
    append!(RC,[(i,) for i in 1:n])
    # add all those edges, which vertices are at most r apart from each other
    for i in 1:n 
        for j in (i+1):n 
            # check the distance of the vertices 
            if norm(collect(S[i]) - collect(S[j])) <= r 
                # add the i,j edge to the complex 
                append!(RC, [(i,j)])
            end
        end
    end
    # add higher dimensional simplices 
    dim = 2
    # simplices of lower dimension
    L  = filter(x -> length(x) == dim, RC)
    while !isempty(L) 
        # increase the dimension by one
        dim += 1
        # a simplex of the dimension dim+1 should be added 
        # if all its faces with dimension dim are contained in RC
        new_simplices = collect(combinations(L, dim))
        new_simplices = [reduce(union,Set(x)) for x in new_simplices]
        filter!(x -> length(x) == dim, new_simplices)
        append!(RC, [tuple(x...) for x in new_simplices])
        # update L 
        L = filter(x -> length(x) == dim, RC)
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