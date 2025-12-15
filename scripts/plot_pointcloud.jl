using CSV
using DataFrames
using Plots

gr()  # backend le plus fiable sous VS Code

function plot_pointcloud(csvfile; title="")
    # 1) Vérif chemin
    if !isfile(csvfile)
        error("File not found: $csvfile (cwd = $(pwd()))")
    end

    # 2) Lecture + aperçu
    df = CSV.read(csvfile, DataFrame)
    println("Loaded $(nrow(df)) points from $csvfile")
    println(first(df, 5))

    # 3) Plot (3D)
    p = scatter(df.x, df.y, df.z;
        markersize=2,
        legend=false,
        title=title,
        aspect_ratio=:equal
    )

    # 4) Affichage VS Code
    display(p)

    # 5) Secours: sauvegarde PNG
    mkpath("results/plots")
    outpng = "results/plots/" * replace(title, " " => "_") * ".png"
    savefig(p, outpng)
    println("Saved PNG: $outpng")

    return p
end

plot_pointcloud("data/torus/samples/torus_0001.csv", title="Torus")
