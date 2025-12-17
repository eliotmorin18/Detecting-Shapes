
using LinearAlgebra

# ============================================================
# 1) Distances
# ============================================================

"""
pairwise_distances(points) -> D

points: N×3 Matrix{Float64}
D[i,j] = ||points[i,:] - points[j,:]||₂
"""
function pairwise_distances(points::Matrix{Float64})
    N = size(points, 1)
    D = zeros(Float64, N, N)

    @inbounds for i in 1:N
        xi = @view points[i, :]
        for j in (i+1):N
            d = norm(xi .- @view(points[j, :]))
            D[i, j] = d
            D[j, i] = d
        end
    end
    return D
end

# ============================================================
# 2) Filtration événementielle (Rips jusqu'à dim=2)
# ============================================================

"""
Un événement de filtration : un simplex apparaît au seuil eps.
- dim = 0 (sommet), 1 (arête), 2 (triangle)
- verts = (i,0,0) / (i,j,0) / (i,j,k), indices triés
"""
struct SimplexEvent
    eps::Float64
    dim::Int
    verts::NTuple{3,Int}
end

vertex_event(i::Int) = SimplexEvent(0.0, 0, (i, 0, 0))
edge_event(i::Int, j::Int, e::Float64) = SimplexEvent(e, 1, (i, j, 0))
tri_event(i::Int, j::Int, k::Int, e::Float64) = SimplexEvent(e, 2, (i, j, k))

"""
build_event_filtration(D; eps_max=2.0, tie_break=true, δ=1e-10) -> events

Construit la liste des événements :
- sommets : eps = 0
- arêtes : eps = D[i,j]
- triangles : eps = max(Dij, Dik, Djk)

eps_max borne la génération pour éviter l'explosion combinatoire.
tie_break ajoute un très petit décalage dépendant de dim pour casser les ex-aequo
(et réduire les cycles H1 de persistance 0).
"""
function build_event_filtration(
    D::AbstractMatrix{<:Real};
    eps_max::Real=2.0,
    tie_break::Bool=true,
    δ::Float64=1e-10
)
    N = size(D, 1)
    events = Vector{SimplexEvent}()
    sizehint!(events, N + div(N*(N-1), 2))

    # Sommets
    @inbounds for i in 1:N
        push!(events, vertex_event(i))
    end

    # Voisinages (pour triangles)
    neighbors = [BitSet() for _ in 1:N]

    # Arêtes
    @inbounds for i in 1:N-1
        for j in i+1:N
            dij = Float64(D[i, j])
            if dij <= eps_max
                push!(events, edge_event(i, j, dij))
                push!(neighbors[i], j)
                push!(neighbors[j], i)
            end
        end
    end

    # Triangles via intersection de voisinages
    @inbounds for i in 1:N-2
        for j in i+1:N-1
            dij = Float64(D[i, j])
            if dij > eps_max
                continue
            end
            common = intersect(neighbors[i], neighbors[j])
            for k in common
                if k > j
                    e = max(dij, Float64(D[i, k]), Float64(D[j, k]))
                    if e <= eps_max
                        push!(events, tri_event(i, j, k, e))
                    end
                end
            end
        end
    end

    # Tie-break (optionnel)
    if tie_break
        @inbounds for t in eachindex(events)
            ev = events[t]
            events[t] = SimplexEvent(ev.eps + ev.dim * δ, ev.dim, ev.verts)
        end
    end

    # Tri final: eps croissant, puis dim croissant
    sort!(events, by = ev -> (ev.eps, ev.dim, ev.verts[1], ev.verts[2], ev.verts[3]))
    return events
end

# ============================================================
# 3) Persistance H0/H1 (réduction colonne mod 2)
# ============================================================

# XOR in-place pour BitSet (différence symétrique)
function xor_bitset!(a::BitSet, b::BitSet)
    @inbounds for x in b
        if x in a
            delete!(a, x)
        else
            push!(a, x)
        end
    end
    return a
end

low(bs::BitSet) = isempty(bs) ? 0 : maximum(bs)

"""
    persistence_pairs(events) -> (pairsH0, pairsH1)

Entrée : events triés (eps croissant, puis dim croissant)
Sortie :
- pairsH0 : Vector{Tuple{Float64,Float64}} (birth, death) avec death=Inf possible
- pairsH1 : Vector{Tuple{Float64,Float64}} (birth, death) avec death=Inf possible
"""
function persistence_pairs(events::Vector{SimplexEvent})
    m = length(events)

    eps_of = Vector{Float64}(undef, m)
    dim_of = Vector{Int}(undef, m)
    verts_of = Vector{NTuple{3,Int}}(undef, m)

    # lookup index dans l'ordre de filtration
    idx0 = Dict{Int,Int}()             # vertex i -> index
    idx1 = Dict{Tuple{Int,Int},Int}()  # edge (i,j) -> index

    for (t, ev) in enumerate(events)
        eps_of[t] = ev.eps
        dim_of[t] = ev.dim
        verts_of[t] = ev.verts

        if ev.dim == 0
            idx0[ev.verts[1]] = t
        elseif ev.dim == 1
            idx1[(ev.verts[1], ev.verts[2])] = t
        end
    end

    # colonnes de bord ∂ (stockées mod 2 dans BitSet)
    cols = Vector{BitSet}(undef, m)
    for t in 1:m
        d = dim_of[t]
        v = verts_of[t]
        if d == 0
            cols[t] = BitSet()
        elseif d == 1
            cols[t] = BitSet((idx0[v[1]], idx0[v[2]]))
        else
            i, j, k = v[1], v[2], v[3]
            cols[t] = BitSet((idx1[(i, j)], idx1[(i, k)], idx1[(j, k)]))
        end
    end

    reduced = Vector{BitSet}(undef, m)
    pivot = Dict{Int,Int}()  # low -> col index

    # pairing info
    death_vertex = Dict{Int,Float64}()  # vertex_index -> death eps (H0)
    birth_edge   = Dict{Int,Float64}()  # edge_index -> birth eps (H1)
    death_edge   = Dict{Int,Float64}()  # edge_index -> death eps (H1)

    # Réduction standard
    for j in 1:m
        col = copy(cols[j])
        lj = low(col)

        while lj != 0 && haskey(pivot, lj)
            xor_bitset!(col, reduced[pivot[lj]])
            lj = low(col)
        end

        reduced[j] = col

        if isempty(col)
            # naissance d'une classe H_dim (ici H1 si dim==1)
            if dim_of[j] == 1
                birth_edge[j] = eps_of[j]
            end
        else
            i = low(col)
            pivot[i] = j

            if dim_of[j] == 1
                # une arête tue une classe H0 (composante)
                death_vertex[i] = eps_of[j]
            elseif dim_of[j] == 2
                # un triangle tue une classe H1 (cycle)
                death_edge[i] = eps_of[j]
            end
        end
    end

    # Construire les paires H0
    pairsH0 = Tuple{Float64,Float64}[]
    for (_vid, vindex) in idx0
        b = eps_of[vindex]
        d = get(death_vertex, vindex, Inf)
        push!(pairsH0, (b, d))
    end

    # Construire les paires H1
    pairsH1 = Tuple{Float64,Float64}[]
    for (eindex, b) in birth_edge
        d = get(death_edge, eindex, Inf)
        push!(pairsH1, (b, d))
    end

    sort!(pairsH0, by = p -> (p[1], p[2]))
    sort!(pairsH1, by = p -> (p[1], p[2]))
    return pairsH0, pairsH1
end

# ============================================================
# 4) Filtres pratiques (garder seulement les grands intervalles)
# ============================================================

"""
    filter_pairs(pairs; min_persistence=0.05, keep_inf=true)

Garde les paires (b,d) telles que:
- d == Inf (si keep_inf)
- ou (d-b) >= min_persistence
"""
function filter_pairs(pairs; min_persistence=0.25, keep_inf=true)
    out = Tuple{Float64,Float64}[]
    for (b, d) in pairs
        if isinf(d)
            keep_inf && push!(out, (b, d))
        else
            (d - b) >= min_persistence && push!(out, (b, d))
        end
    end
    return out
end
