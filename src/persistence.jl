# src/persistence.jl
module Persistence

export create_custom_filtration

using LinearAlgebra, Plots, Ripserer, PersistenceDiagrams 

include("filtration.jl")
import .Filtration

""" Computing persistence homology wit hthe help of useful packages """

function create_custom_filtration(filtration)
    # from filtration to the custom filtration that can be used as input in ripserer 

    # INPUT 
    # filtration 

    # OUTPUT 
    # custom 

    # K = Custom([(1,) => 0, (3,) => 0, (2,4) => 1, (1,4) => 2, (2,3) => 2, (3,4) => 3, (1,2) => 4, (1,2,4) => 5])
    k = sort(collect(keys(filtration)))
    flt = Dict( findfirst(x -> x==r, k)  => [sort(sx) for sx in filtration[r]] for r in keys(filtration))
    c_pair = [sx => 1 for sx in flt[1]]
    for i in 2:length(k) 
        add_simplices = setdiff(flt[i], flt[i-1])
        add_pair = [sx => i for sx in add_simplices]
        append!(c_pair, add_pair)
    end 

    return Custom(c_pair)
end 

function create_wasserstein_matrix(pointclouds, filtration_type, print_output)
    # crate a matrix containing the wasserstein distances for all the pointclouds

    # INPUT 
    # pointclouds = Vecor, each entry is a pointcloud 
    # filtration_type = "Rips", "Cech", "Ripserer" 

    # ASSUMPTIONS 
    # all pointclouds have the same dimension

    # OUTPUT
    # Ds = vector containing the distance matrices of persistence diagrams for H_i i in 0:(dim-1) 
    # Ds[hi] = distance matrix for hi in 1:dim == 0:(dim-1)
    # Ds[hi][i,j] = Ds[hi][j,i] = Wasserstein Distance between persistence diagrams of the point clouds

    if print_output
        println("Create Wasserstein Matrix with the filtration_type ", filtration_type)
    end 

    n_pointclouds = length(pointclouds)
    dim = length(pointclouds[1][1])

    if print_output
        println("n_pointclouds = ", n_pointclouds)
        println("dim = ", dim)
    end 

    results = [] 
    Ds = [zeros(n_pointclouds, n_pointclouds) for i in 1:(dim-1)]

    for i in 1:n_pointclouds 
        # for each pointcloud 
        # 1) if filtration_type != "Ripserer" : compute filtration, convert the filtration to a custom filtration 
        # 2) compute homology using ripserer 
        # 3) get Wasserstein distance to all the other pointclouds 
        # 4) save distances in the matrix 
        S = pointclouds[i]

        if print_output
            println("Compute homology of the pointcloud ", i)
        end 

        if filtration_type == "Ripserer"
            result = ripserer(S; dim_max = dim-1)
        elseif filtration_type == "Rips"
            # different radii to use 
            R = Filtration.quantile_distances(S)
            
            if print_output 
                println("start to compute the filtration")
            end 
            
            filtration = Filtration.rips_filtration(S, R)
            
            if print_output
                println("filtration computed")
            end 

            custom_filtration = create_custom_filtration(filtration)
            result = ripserer(custom_filtration; alg = :homology)
        elseif filtration_type == "Cech"
            # different radii to use 
            R = Filtration.quantile_distances(S)
            R = 1/2 .* R

            if print_output 
                println("start to compute the filtration")
            end 
            
            filtration = Filtration.cech_filtration(S, R)
            
            if print_output
                println("filtration computed")
            end 

            custom_filtration = create_custom_filtration(filtration)
            result = ripserer(custom_filtration; alg = :homology)
        else 
            println("ERROR filtration type not equal to one of the following: Ripserer, Cech, Rips")
            return nothing 
        end 

        if print_output
            println("result computed")
        end 
        
        # save the result 
        push!(results, result)
        
        # compute and save the distances
        if print_output
            println("start to compute the distances")
        end

        for j in 1:(i-1)
            for hi in 1:(dim-1) 
                d1 = Wasserstein()(result[hi], results[j][hi])
                d2 = Wasserstein()(results[j][hi], result[hi])
                # the computed distances are not exactly the same, if e.g. one does not have any entries in H_2 while the second one has, the distance is zero
                # interchanging the arguments give a positive distance 
                # Thus we take the maximum 
                d = maximum([d1,d2])
                Ds[hi][i,j] = d
                Ds[hi][j,i] = d 
            end 
        end
        if print_output
            println("distances computed")
        end 
    end 

    

    

    return Ds
end 






end # end module 