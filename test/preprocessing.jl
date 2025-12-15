# test/preprocessing.jl

using Statistics
using LinearAlgebra

# Include the preprocessing module
include("../src/preprocessing.jl")
using .Preprocessing

# Path to your test CSV
path = "data/circle/samples/circle_0001.csv"

println("=== PREPROCESSING TEST ===")
println("File: ", path)

# Run preprocessing
P = prepare_points(path)

println("\nNumber of points: ", size(P, 1))

# ---- Test 1: mean should be ~ 0 ----
mean_point = vec(mean(P, dims=1))
println("\nMean of the point cloud (should be close to 0):")
println(mean_point)

# ---- Test 2: farthest point should be at distance 1 ----
norms = sqrt.(sum(P.^2, dims=2))[:]
rmax = maximum(norms)
idx_max = argmax(norms)

println("\nMaximum distance to origin (should be 1):")
println(rmax)

println("\nIndex of farthest point:")
println(idx_max)

println("\nCoordinates of the farthest point:")
println(P[idx_max, :])

# ---- Show a few sample points ----
println("\nFirst 5 points after preprocessing:")
println(P[1:min(5, size(P,1)), :])

println("\n=== TEST FINISHED ===")
