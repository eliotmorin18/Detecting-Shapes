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

# Fonction de perturbation directionnelle f(u)
# u = (ux,uy,uz) unitaire, sortie ~ [-1,1] (en pratique)
function direction_perturbation(u::AbstractVector{<:Real}; k1=3.0, k2=5.0, k3=7.0)
    ux, uy, uz = u
    return 0.5*sin(k1*ux) + 0.35*cos(k2*uy) + 0.25*sin(k3*uz)
end

function generate_potato_pointcloud(rng::AbstractRNG; N::Int=2000,
                                    R::Float64=3.0, α::Float64=0.25,
                                    σ::Float64=0.01, translate_scale::Float64=5.0)
    X = zeros(Float64, N, 3)
    noise = Normal(0, σ)

    for i in 1:N
        # direction uniforme sur la sphère
        u = randn(rng, 3)
        u /= norm(u)

        f = direction_perturbation(u)
        rmax = R * (1 + α * f)
        # sécurité (au cas où rmax devienne trop petit)
        rmax = max(rmax, 0.2R)

        # rayon uniforme en volume : r = rmax * cbrt(U)
        r = rmax * cbrt(rand(rng))

        p = r .* u
        p[1] += rand(rng, noise)
        p[2] += rand(rng, noise)
        p[3] += rand(rng, noise)

        X[i, :] .= p
    end

    # rotation + translation
    Rot = random_rotation(rng)
    X = X * Rot'

    d = (rand(rng, 3) .- 0.5) .* (2 * translate_scale)
    X .+= d'

    return X, Rot, d
end

function make_potato_dataset(; K::Int=30, N::Int=2000,
                             R::Float64=3.0, α::Float64=0.25,
                             σ::Float64=0.01, translate_scale::Float64=5.0,
                             outdir::String="data/potato")
    samples_dir = joinpath(outdir, "samples")
    mkpath(samples_dir)

    meta = DataFrame(
        sample = String[],
        shape = String[],
        N = Int[],
        R = Float64[],
        alpha = Float64[],
        sigma = Float64[],
        translate_scale = Float64[],
        seed = Int[]
    )

    for i in 1:K
        seed = 60_000 + i
        rng = MersenneTwister(seed)

        X, _, _ = generate_potato_pointcloud(rng; N=N, R=R, α=α, σ=σ,
                                             translate_scale=translate_scale)

        filename = @sprintf("potato_%04d.csv", i)
        filepath = joinpath(samples_dir, filename)

        df = DataFrame(x=X[:,1], y=X[:,2], z=X[:,3])
        CSV.write(filepath, df)

        push!(meta, (filename, "potato", N, R, α, σ, translate_scale, seed))
    end

    CSV.write(joinpath(outdir, "metadata.csv"), meta)
    println("✅ Potato (perturbed 3-disc) dataset created in: ", outdir)
    println("   - samples: ", K, " CSV files")
    println("   - metadata.csv written")
end

# ---- Lance la génération ----
make_potato_dataset(K=30, N=400, R=3.0, α=0.25, σ=0.01, translate_scale=5.0)
