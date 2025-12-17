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
    
    for i in 1:n 
        for j in (i+1):n 
            # check if the balls of radius r would intersect, i.e. dist <= 2r
            if norm(collect(S[i]) - collect(S[j])) <= 2r 
                # add the i,j edge to the complex 
                append!(CC, [(i,j)])
            end
        end
    end
    # add higher dimensional simplices 
    dim = 2 
    L  = filter(x -> length(x) == dim, CC)
    while !isempty(L) 
        dim += 1 
        # possible simplices to add are only those, where all the faces are contained in CC
        new_simplices = collect(combinations(L, dim))
        new_simplices = [reduce(union,Set(x)) for x in new_simplices]
        filter!(x -> length(x) == dim, new_simplices)
        
        # for those simplices check if balls intersect 
        filter!(x -> minimal_enclosing_ball([collect(p) for p in S[x]])[1] <= r, new_simplices)

        append!(CC, [tuple(x...) for x in new_simplices])
        L  = filter(x -> length(x) == dim, CC)
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
        origin = S[1] + 1/2* (S[2] .- S[1])
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
            A[i-1, :] = 2 * (S[i] .- S[1])
            b[i-1] = norm(S[i])^2 - norm(S[1])^2
        end
        origin = A\b 
        return norm(origin - S[1]), origin
    end 
end 