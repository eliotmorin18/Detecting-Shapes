using Random
using Distributions
using LinearAlgebra
using CSV
using DataFrames
using Printf

# --- rotation aléatoire (matrice orthogonale det=+1) ---
function random_rotation(rng::AbstractRNG)
    M = randn(rng, 3, 3)
    Q, _ = qr(M)
    Q = Matrix(Q)
    if det(Q) < 0
        Q[:, 1] .*= -1
    end
    return Q
end

function generate_line_segment_pointcloud(rng::AbstractRNG; N::Int=300, L::Float64=10.0, σ::Float64=0.05,
                                         translate_scale::Float64=5.0)
    A = [-L/2, 0.0, 0.0]
    B = [ L/2, 0.0, 0.0]

    t = rand(rng, N)  # Uniform(0,1)
    X = zeros(Float64, N, 3)

    noise = Normal(0, σ)

    for i in 1:N
        P = A .+ t[i] .* (B .- A)
        # bruit "thin" (épaisseur) en y,z
        P[2] += rand(rng, noise)
        P[3] += rand(rng, noise)
        X[i, :] .= P
    end

    # rotation + translation
    R = random_rotation(rng)
    X = X * R'  # rotation

    d = (rand(rng, 3) .- 0.5) .* (2 * translate_scale)  # dans [-scale, scale]
    X .+= d'  # translation

    return X, R, d
end

function make_line_segment_dataset(; K::Int=50, N::Int=300, L::Float64=10.0, σ::Float64=0.05,
                                   translate_scale::Float64=5.0, outdir::String="data/line_segment")
    samples_dir = joinpath(outdir, "samples")
    mkpath(samples_dir)

    meta = DataFrame(
        sample = String[],
        shape = String[],
        N = Int[],
        L = Float64[],
        sigma = Float64[],
        translate_scale = Float64[],
        seed = Int[]
    )

    for i in 1:K
        seed = 10_000 + i
        rng = MersenneTwister(seed)

        X, _, _ = generate_line_segment_pointcloud(rng; N=N, L=L, σ=σ, translate_scale=translate_scale)

        filename = @sprintf("line_%04d.csv", i)
        filepath = joinpath(samples_dir, filename)

        df = DataFrame(x=X[:,1], y=X[:,2], z=X[:,3])
        CSV.write(filepath, df)

        push!(meta, (filename, "line_segment", N, L, σ, translate_scale, seed))
    end

    CSV.write(joinpath(outdir, "metadata.csv"), meta)
    println("Dataset created in: ", outdir)
    println("   - samples: ", K, " CSV files")
    println("   - metadata.csv written")
end

# ---- Lance la génération ----
make_line_segment_dataset(K=50, N=300, L=10.0, σ=0.05, translate_scale=5.0)
