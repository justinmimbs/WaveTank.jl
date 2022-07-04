module MonaiValley

using CSV
using Interpolations: interpolate, scale, BSpline, Cubic, Flat, OnCell
using GLMakie
using Neowave: g, Grid, Model, run!, Results

import ..to_path
import ..gridindex, ..sample
import ..plot_scene, ..header!

function model()
    csv = CSV.File(to_path("in/monaivalley_bathymetry.csv"); types=Float64)
    dx = csv[2].y - csv[1].y
    rx = csv[1].x:dx:csv[end].x
    ry = csv[1].y:dx:csv[end].y
    z = reshape(csv.z, length(ry), length(rx))
    #
    h = 0.5dx
    nx = length(rx)
    ny = length(ry)
    grid = Grid((rx[1] - h, rx[end] + h), (ry[1] - h, ry[end] + h), (nx, ny))
    function bathymetry(x, y)
        i, j = round.(Int, (y, x) ./ dx) .+ 1
        z[i, j]
    end
    Model(;
        grid,
        h = bathymetry,
        bcx = :wall,
        bcy = :wall,
        dt = 0.005,
        n = 0.012,
    )
end

function inputwave(; dt)
    csv = CSV.File(to_path("in/monaivalley_inputwave.csv"); types=Float64)
    s, eta = csv.s, csv.eta
    #
    s1 = range(s[1], s[end], length=length(s))
    etaf = scale(interpolate(eta, BSpline(Cubic(Flat(OnCell())))), s1)
    s2 = range(s[1], s[end], step=dt)
    [ etaf(s) for s in s2 ]
end

function run(name="monaivalley")
    outfile = to_path("out/$name.jld2")
    if isfile(outfile)
        @info "file exists: $(outfile)"
    else
        m = model()
        waveinput = inputwave(; m.dt)
        run!(m, outfile; seconds=27.0, frequency=20, output=(; eta=:eta), waveinput)
        @info "file created: $(outfile)"
    end
end

#

gauge_positions() = [ (4.521, 1.196), (4.521, 1.696), (4.521, 2.196) ]

function plot_setup()
    (; grid, h) = model()
    gauges = gauge_positions()
    dt = 0.1
    wave = inputwave(; dt)

    fig = Figure()
    header!(fig[1, 1], "Monai Valley, setup")

    # bathymetry
    ax = Axis(fig[2, 1]; title="bathymetry", xlabel="x (m)", ylabel="y (m)",
        aspect=DataAspect(),
        limits=(grid.xf[1], grid.xf[end], grid.yf[1], grid.yf[end]),
        xticks=0:5,
        yticks=0:3,
    )
    # terrain contours
    l, u = extrema(h)
    dz = 0.01 # (u - l) / 27
    levels = vcat(-dz:-dz:l, dz:dz:u)
    contour!(ax, grid.xc, grid.yc, h; levels, color=:gray)
    # nominal waterline
    contour!(ax, grid.xc, grid.yc, h; levels=[0.0], color=:black, linewidth=1.0)
    # gauges
    scatter!(ax, gauges; color=:black, marker=:circle, markersize=12, label="gauges")
    text!(ax, gauges; text=string.(1:length(gauges)), align=(:right, :center), offset=(-6, 0))
    axislegend(ax; position=(:right, :bottom))

    # wave input
    axw = Axis(fig[3, 1]; title="wave input", xlabel="time (s)", ylabel="surface elevation (cm)",
        limits=(0, 25, -2, 2),
    )
    # hidespines!(axw)
    s = 0.0:dt:((length(wave) - 1) * dt)
    lines!(axw, s, wave .* 100; color=:black)
    rowsize!(fig.layout, 3, Auto(0.3))

    fig
end

function plot_timeseries(name="monaivalley")
    positions = gauge_positions()
    n = length(positions)
    # simulation data
    filename = to_path("out/$name.jld2")
    (; grid) = first(Results(filename))
    eta_sim = sample(filename, gridindex.(Ref(grid), positions))
    t = length(eta_sim[1])
    # lab data
    csv = CSV.File(to_path("in/monaivalley_wavegages.csv"); types=Float64)
    seconds = csv.s[1:t]
    eta_lab = [ c.column[1:t] for c in csv.columns[1 .+ (1:n)] ]

    fig = Figure()
    header!(fig[1, 1], "Monai Valley, timeseries")

    limits = (seconds[1], seconds[end], -1.2, 5.2)
    for i in 1:n
        lab, sim = eta_lab[i], eta_sim[i]
        ax = Axis(fig[i + 1, 1]; title="gauge $i", limits, xlabel="time (s)", ylabel="eta (cm)")
        # hidespines!(ax)
        lines!(ax, seconds, lab; color=:black, label="laboratory")
        lines!(ax, seconds, sim .* 100; color=:dodgerblue, label="model", linewidth=3)
        if i == 1; axislegend(ax) end
    end

    fig
end

function plot_results(name="monaivalley")
    filename = to_path("out/$name.jld2")
    plot_scene(Results(filename), "Monai Valley"; zscale=5.0, dz=0.01)
end

end # module
