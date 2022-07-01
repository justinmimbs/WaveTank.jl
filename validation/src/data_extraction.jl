using JLD2: jldopen

function lineindex(xs, x)
    a = (x - xs[1]) / (xs[end] - xs[1])
    round(Int, a * (length(xs) - 1)) + 1
end

function gridindex(grid, (x, y))
    lineindex(grid.xc, x), lineindex(grid.yc, y)
end

# collect timeseries data at one grid point
function sample(filepath, (i, j))
    jldopen(filepath, "r") do file
        ts = sort!(parse.(Int, keys(file["output/eta"])))
        [ file["output/eta/$t"][i, j] for t in ts ]
    end
end

# collect timeseries data at multiple grid points
function sample(filepath, indicies::Vector)
    jldopen(filepath, "r") do file
        ts = sort!(parse.(Int, keys(file["output/eta"])))
        data = [ zeros(length(ts)) for _ in indicies ]
        for (i, t) in enumerate(ts)
            eta = file["output/eta/$t"]
            for (series, ij) in zip(data, indicies)
                series[i] = eta[ij...]
            end
        end
        data
    end
end
