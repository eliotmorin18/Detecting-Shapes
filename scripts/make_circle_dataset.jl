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

function generate_circle_pointcloud(rng::AbstractRNG; N::Int=400, R::Float64=3.0,
                                    σ::Float64=0.03, translate_scale::Float64=5.0)
    # Cercle de base dans le plan z=0
    θ = 2π .* rand(rng, N)
    X = zeros(Float64, N, 3)

    noise = Normal(0, σ)

    for i in 1:N
        x = R * cos(θ[i])
        y = R * sin(θ[i])
        z = 0.0

        # bruit (petit) en 3D
        x += rand(rng, noise)
        y += rand(rng, noise)
        z += rand(rng, noise)

        X[i, :] .= (x, y, z)
    end

    # rotation + translation
    Rot = random_rotation(rng)
    X = X * Rot'

    d = (rand(rng, 3) .- 0.5) .* (2 * translate_scale)
    X .+= d'

    return X, Rot, d
end

function make_circle_dataset(; K::Int=50, N::Int=400, R::Float64=3.0, σ::Float64=0.03,
                             translate_scale::Float64=5.0, outdir::String="data/circle")
    samples_dir = joinpath(outdir, "samples")
    mkpath(samples_dir)

    meta = DataFrame(
        sample = String[],
        shape = String[],
        N = Int[],
        R = Float64[],
        sigma = Float64[],
        translate_scale = Float64[],
        seed = Int[]
    )

    for i in 1:K
        seed = 20_000 + i
        rng = MersenneTwister(seed)

        X, _, _ = generate_circle_pointcloud(rng; N=N, R=R, σ=σ, translate_scale=translate_scale)

        filename = @sprintf("circle_%04d.csv", i)
        filepath = joinpath(samples_dir, filename)

        df = DataFrame(x=X[:,1], y=X[:,2], z=X[:,3])
        CSV.write(filepath, df)

        push!(meta, (filename, "circle", N, R, σ, translate_scale, seed))
    end

    CSV.write(joinpath(outdir, "metadata.csv"), meta)
    println("✅ Circle dataset created in: ", outdir)
    println("   - samples: ", K, " CSV files")
    println("   - metadata.csv written")
end

# ---- Lance la génération ----
make_circle_dataset(K=50, N=80, R=3.0, σ=0.03, translate_scale=5.0)
