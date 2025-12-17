using LinearAlgebra
using Printf

# --- Load project code
include(joinpath(@__DIR__, "..", "src", "preprocessing.jl"))  # defines prepare_points(...)
include(joinpath(@__DIR__, "..", "src", "filtration.jl"))     # defines pairwise_distances, build_event_filtration, persistence_pairs, filter_pairs

# --- Parameters
csv_path = joinpath(@__DIR__, "..", "data", "potato", "samples", "potato_0004.csv")
NKEEP    = 80
EPS_MAX  = 1.0        # start small to avoid combinatorial explosion

# --- Pipeline
P_all = prepare_points(csv_path)
n_all = size(P_all, 1)
println("CSV: ", csv_path)
println("Total points loaded: ", n_all)

n_keep = min(NKEEP, n_all)
P = P_all[1:n_keep, :]
println("Points kept: ", size(P, 1), " (first ", n_keep, ")")

D = pairwise_distances(P)

events = build_event_filtration(D; eps_max=EPS_MAX, tie_break=true)  # tie_break helps reduce zero-length H1
println("Events built: ", length(events), "  (eps_max=", EPS_MAX, ")")

# --- Pretty-print events
function print_event(ev)
    if ev.dim == 0
        @printf("eps = %.6f | vertex   %d\n", ev.eps, ev.verts[1])
    elseif ev.dim == 1
        @printf("eps = %.6f | edge     (%d,%d)\n", ev.eps, ev.verts[1], ev.verts[2])
    else
        @printf("eps = %.6f | triangle (%d,%d,%d)\n", ev.eps, ev.verts[1], ev.verts[2], ev.verts[3])
    end
end

println("\n==== First 20 events ====")
for ev in events[1:min(20, length(events))]
    print_event(ev)
end

println("\n==== Last 20 events ====")
for ev in events[max(1, length(events)-19):end]
    print_event(ev)
end

# --- Persistence
pairsH0, pairsH1 = persistence_pairs(events)
println("\nPersistence summary")
println("  |H0| = ", length(pairsH0), "   (expect ~", n_keep, " intervals; 1 should be Inf)")
println("  |H1| = ", length(pairsH1))

# sanity checks (prints, not tests)
n_inf_H0 = count(p -> isinf(p[2]), pairsH0)
println("  H0 intervals with death=Inf: ", n_inf_H0)

# show some H1 intervals (filtered)
pairsH1_big = filter_pairs(pairsH1; min_persistence=0.2, keep_inf=true)
println("  H1 intervals with persistence >= 0.45 (or Inf): ", length(pairsH1_big))

println("\n==== First 10 H1 (filtered) ====")
for p in pairsH1_big[1:min(10, length(pairsH1_big))]
    @printf("(birth=%.6f, death=%s, pers=%s)\n",
            p[1],
            isinf(p[2]) ? "Inf" : @sprintf("%.6f", p[2]),
            isinf(p[2]) ? "Inf" : @sprintf("%.6f", p[2]-p[1]))
end
