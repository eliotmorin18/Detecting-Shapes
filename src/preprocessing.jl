
export load_pointcloud_csv, center_points, normalize_scale, prepare_points

using CSV
using DataFrames
using Statistics
using LinearAlgebra

"""
    load_pointcloud_csv(path::AbstractString) -> Matrix{Float64}

Input:
- path: path to a CSV file containing 3 columns (x,y,z), with or without a header.

Output:
- P: N×3 Matrix{Float64} where each row is a point [x y z].

Loads a point cloud from CSV and returns a numeric matrix ready for processing.
"""


function load_pointcloud_csv(path::AbstractString)::Matrix{Float64}
    df = CSV.read(path, DataFrame)

    # If the file has more than 3 columns, we keep only the first 3.
    # If it has exactly 3, we use all of them.
    if ncol(df) < 3
        error("CSV must contain at least 3 columns for x,y,z.")
    end

    df3 = df[:, 1:3]
    P = Matrix{Float64}(df3)
    return P
end


"""
    center_points(P::AbstractMatrix{<:Real}) -> Matrix{Float64}

Input:
- P: N×3 matrix (rows are points)

Output:
- Pc: N×3 centered matrix (mean is moved to the origin)

Centers the cloud by subtracting the mean of each coordinate.
"""
function center_points(P::AbstractMatrix{<:Real})::Matrix{Float64}
    Pf = Matrix{Float64}(P)
    c = vec(mean(Pf, dims=1))              # center (x̄, ȳ, z̄)
    Pc = Pf .- reshape(c, 1, :)            # subtract from every row
    return Pc
end


"""
    normalize_scale(P::AbstractMatrix{<:Real}; target=1.0, eps=1e-12) -> Matrix{Float64}

Input:
- P: N×3 matrix (usually already centered)
- target: desired maximum radius after scaling (default: 1.0)
- eps: safety threshold to avoid division by ~0

Output:
- Pn: N×3 normalized matrix

Scales the cloud so that max(||p||) == target.
This makes ε comparable across different objects.
"""
function normalize_scale(
    P::AbstractMatrix{<:Real};
    target::Real = 1.0,
    eps::Real = 1e-12
)::Matrix{Float64}
    Pf = Matrix{Float64}(P)
    norms = sqrt.(sum(Pf.^2, dims=2))[:]   # Euclidean norm per point
    rmax = maximum(norms)

    if rmax < eps
        # Degenerate case: all points are (almost) identical
        return copy(Pf)
    end

    Pn = Pf .* (target / rmax)
    return Pn
end


"""
    prepare_points(path::AbstractString; do_center=true, do_normalize=true) -> Matrix{Float64}

Input:
- path: CSV file path
- do_center: apply centering (default true)
- do_normalize: apply normalization (default true)

Output:
- Pout: processed N×3 Matrix{Float64}

Convenience function: load -> center -> normalize.
"""
function prepare_points(
    path::AbstractString;
    do_center::Bool = true,
    do_normalize::Bool = true
)::Matrix{Float64}
    P = load_pointcloud_csv(path)
    if do_center
        P = center_points(P)
    end
    if do_normalize
        P = normalize_scale(P)
    end
    return P
end

