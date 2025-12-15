using CSV
using DataFrames
using Plots
using Printf

gr()  # backend fiable

# --- Plot + save (PNG) ---
function save_pointcloud_plot(csvfile::String; outpng::String, title::String="")
    df = CSV.read(csvfile, DataFrame)

    p = scatter(
        df.x, df.y, df.z;
        markersize=2,
        legend=false,
        title=title,
        aspect_ratio=:equal
    )

    mkpath(dirname(outpng))
    savefig(p, outpng)
    println("✅ Saved: $outpng")
end

# --- Utilitaire: prend les k premiers CSV d'un dossier samples/ ---
function first_k_csv(samples_dir::String, k::Int)
    files = filter(f -> endswith(lowercase(f), ".csv"), readdir(samples_dir; join=true))
    sort!(files)
    return files[1:min(k, length(files))]
end

# --- Paramètres ---
DATA_DIR = "data"
OUT_DIR  = "results/plots"
K_PER_SHAPE = 2   # mets 1 si tu veux un seul exemple par objet

mkpath(OUT_DIR)

# Parcourt data/<shape>/samples/*.csv
for shape in sort(readdir(DATA_DIR))
    samples_dir = joinpath(DATA_DIR, shape, "samples")
    if !isdir(samples_dir)
        continue
    end

    csvs = first_k_csv(samples_dir, K_PER_SHAPE)
    if isempty(csvs)
        @warn "No CSV found in $samples_dir"
        continue
    end

    for (i, csvfile) in enumerate(csvs)
        title = "$(shape) $(i)"
        outpng = joinpath(OUT_DIR, @sprintf("%s_%02d.png", shape, i))
        save_pointcloud_plot(csvfile; outpng=outpng, title=title)
    end
end

println("🎉 Done. Check: $OUT_DIR")
