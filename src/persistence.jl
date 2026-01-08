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

function create_persistence_diagram(pointclouds, filtration_type, print_output)
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
            result = ripserer(custom_filtration; alg = :homology, dim_max = dim-1)
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
    end 
    return results
end 

function create_wasserstein_matrix(diagrams)
    # crate a matrix containing the wasserstein distances for all the pointclouds
    n_dim = length(diagrams[1])
    n_diagrams = length(diagrams)
    Ds = [zeros(n_diagrams, n_diagrams) for i in 1:n_dim]

    for (i,diagram) in enumerate(diagrams)
        for j in 1:(i-1)
            for hi in 1:n_dim
                d1 = Wasserstein()(diagram[hi], diagrams[j][hi])
                d2 = Wasserstein()(diagrams[j][hi], diagram[hi])
                # the computed distances are not exactly the same, if e.g. one does not have any entries in H_2 while the second one has, the distance is zero
                # interchanging the arguments give a positive distance 
                # Thus we take the maximum 
                d = maximum([d1,d2])
                Ds[hi][i,j] = d
                Ds[hi][j,i] = d 
            end 
        end
    end 
    return Ds
end 

function persistence_silhoutte(diagram; T = 100, p=1)
    # computes the persistence silhouette using the definitions from this paper https://arxiv.org/pdf/1312.0308

    # INPUT 
    # diagram = persistence diagram, out put from ripserer 
    # T = 100 integer, indicating the number of timepoints on time axis 
    # p = 1 or another number, for weight computation 

    # OUTPUT 
    # silhouette = vector representing the value of the silhuette function along the time 

    if length(diagram) == 0 
        return zeros(T)
    end 

    # the bars of the barcode with all the entries which dies in finite time
    bars = [(d.birth, d.death) for d in diagram if d.death < Inf]

    # the weights for the persistence silhouette are given by (death - birth)^p 
    weights = [(bar[2] - bar[1])^p for bar in bars]

    # the time vector from first birth till last finite death
    ts = range(minimum(bar[1] for bar in bars), maximum(bar[2] for bar in bars), length = T)

    # initiate silhouette vector
    sil = zeros(length(ts))

    for ((b, d), w) in zip(bars, weights)
        # b = birth, d = death, w = weight 
        mid = (b + d)/2
        for (j,t) in enumerate(ts)
            if b <= t <= mid 
                sil[j] += w* (t-b)
            elseif mid < t <= d
                sil[j] += w* (d - t)
            end 
        end 
    end 
    
    sil = sil ./ sum(weights) 
    return sil
end 

function create_persistence_silhoutte(diagrams, weights)
    # creates the concatenated silhouettes for all the diagrams 

    # INPUT 
    # diagrams = vector, each entry is a list of the persistence diagrams for that pointcloud in the different dimensions (for H_0, H_1...)

    # OUTPUT 
    # silhouettes =  list contatining all the concatenated silhouttes for all the pointclouds 

    # concatenate all the persistence silhuettes
    silhouettes = []
    for all_diag in diagrams
        # concatenated silhouettes 

        sil_total = []
        for (i,diag) in enumerate(all_diag) 

            # normalize in each dimension separatly 
            sil = persistence_silhoutte(diag) 
            if !(norm(sil) == 0)
                sil ./ norm(sil)
            end 
            # concatenate 
            append!(sil_total,weights[i] .* sil)
        end 
        # add concatenated sequence to all the others 
        push!(silhouettes, sil_total)
    end 
    return silhouettes
end 
end # end module 