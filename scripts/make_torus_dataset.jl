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

function generate_torus_pointcloud(rng::AbstractRNG; N::Int=1000,
                                   R::Float64=4.0, r::Float64=1.2,
                                   σ::Float64=0.02,
                                   translate_scale::Float64=5.0)
    X = zeros(Float64, N, 3)
    noise = Normal(0, σ)

    θ = 2π .* rand(rng, N)
    φ = 2π .* rand(rng, N)

    for i in 1:N
        x = (R + r*cos(θ[i])) * cos(φ[i])
        y = (R + r*cos(θ[i])) * sin(φ[i])
        z = r * sin(θ[i])

        # bruit
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

function make_torus_dataset(; K::Int=40, N::Int=1000,
                            R::Float64=4.0, r::Float64=1.2,
                            σ::Float64=0.02,
                            translate_scale::Float64=5.0,
                            outdir::String="data/torus")
    samples_dir = joinpath(outdir, "samples")
    mkpath(samples_dir)

    meta = DataFrame(
        sample = String[],
        shape = String[],
        N = Int[],
        R = Float64[],
        r = Float64[],
        sigma = Float64[],
        translate_scale = Float64[],
        seed = Int[]
    )

    for i in 1:K
        seed = 70_000 + i
        rng = MersenneTwister(seed)

        X, _, _ = generate_torus_pointcloud(rng; N=N, R=R, r=r,
                                            σ=σ,
                                            translate_scale=translate_scale)

        filename = @sprintf("torus_%04d.csv", i)
        filepath = joinpath(samples_dir, filename)

        df = DataFrame(x=X[:,1], y=X[:,2], z=X[:,3])
        CSV.write(filepath, df)

        push!(meta, (filename, "torus", N, R, r, σ, translate_scale, seed))
    end

    CSV.write(joinpath(outdir, "metadata.csv"), meta)
    println("✅ Torus dataset created in: ", outdir)
    println("   - samples: ", K, " CSV files")
    println("   - metadata.csv written")
end

# ---- Lance la génération ----
make_torus_dataset(K=40, N=300, R=4.0, r=1.2, σ=0.02, translate_scale=5.0)
