# src/visualization.jl
module Visualization

export pointcloud_3d

using Plots 

function pointcloud_3d(S; check_npoints = true)
    # plots a pointcloud in 3d 

    # INPUT 
    # S = point cloud 

    # EXAMPLE INPUT 
    # S = [(1,2,3), (4,2,1), (2,2,2)]

    # OUTPUT 
    # plot 
    markersize = 4
    n = length(S)
    if (n > 100) & check_npoints 
        println("There are a lot of points for the plot.
         If you really want to make the plot, set the argument check_npoints = false")
        return nothing 
    end 

    x = first.(S)
    y = getindex.(S, 2)
    z = last.(S)

    pl = scatter(
        x, y, z;
        seriestype = :scatter,
        markersize = markersize,
        xlabel = "x",
        ylabel = "y",
        zlabel = "z",
        legend = false
    )

end 
end # end module 