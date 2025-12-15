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

# Échantillonnage uniforme dans un disque de rayon R
# Astuce: r = R*sqrt(u) pour être uniforme en aire
function generate_disk_pointcloud(rng::AbstractRNG; N::Int=600, R::Float64=3.0,
                                  σ_xy::Float64=0.03, σ_z::Float64=0.01,
                                  translate_scale::Float64=5.0)
    u = rand(rng, N)          # pour le rayon
    θ = 2π .* rand(rng, N)    # pour l'angle
    r = R .* sqrt.(u)

    X = zeros(Float64, N, 3)

    noise_xy = Normal(0, σ_xy)
    noise_z  = Normal(0, σ_z)

    for i in 1:N
        x = r[i] * cos(θ[i]) + rand(rng, noise_xy)
        y = r[i] * sin(θ[i]) + rand(rng, noise_xy)
        z = 0.0 + rand(rng, noise_z)   # disque "mince"
        X[i, :] .= (x, y, z)
    end

    # rotation + translation
    Rot = random_rotation(rng)
    X = X * Rot'

    d = (rand(rng, 3) .- 0.5) .* (2 * translate_scale)
    X .+= d'

    return X, Rot, d
end

function make_disk_dataset(; K::Int=50, N::Int=600, R::Float64=3.0,
                           σ_xy::Float64=0.03, σ_z::Float64=0.01,
                           translate_scale::Float64=5.0,
                           outdir::String="data/disk")
    samples_dir = joinpath(outdir, "samples")
    mkpath(samples_dir)

    meta = DataFrame(
        sample = String[],
        shape = String[],
        N = Int[],
        R = Float64[],
        sigma_xy = Float64[],
        sigma_z = Float64[],
        translate_scale = Float64[],
        seed = Int[]
    )

    for i in 1:K
        seed = 30_000 + i
        rng = MersenneTwister(seed)

        X, _, _ = generate_disk_pointcloud(rng; N=N, R=R, σ_xy=σ_xy, σ_z=σ_z,
                                           translate_scale=translate_scale)

        filename = @sprintf("disk_%04d.csv", i)
        filepath = joinpath(samples_dir, filename)

        df = DataFrame(x=X[:,1], y=X[:,2], z=X[:,3])
        CSV.write(filepath, df)

        push!(meta, (filename, "disk", N, R, σ_xy, σ_z, translate_scale, seed))
    end

    CSV.write(joinpath(outdir, "metadata.csv"), meta)
    println("✅ Disk dataset created in: ", outdir)
    println("   - samples: ", K, " CSV files")
    println("   - metadata.csv written")
end

# ---- Lance la génération ----
make_disk_dataset(K=50, N=600, R=3.0, σ_xy=0.03, σ_z=0.01, translate_scale=5.0)
