module ConicalIsland

using CSV
using GLMakie
using Printf: @sprintf
using Neowave: g, Grid, Model, run!, Results, load, foldtime

import ..to_path
import ..sech2wave, ..particle_velocity, ..truncatedcone
import ..gridindex, ..sample
import ..plot_scene, ..header!

function model(res=1)
    h = 0.32
    a = 0.181 * h
    k = sqrt(3 * a / (4 * h^3))
    wave = sech2wave(; a, x0=7, k)
    speed = particle_velocity(wave, a, h)
    basin = truncatedcone(13, 15, 3.6, 1.1, 0.625)
    Model(
        grid=Grid(0:25, 0:30, (500 ÷ res, 600 ÷ res)),
        h=(x, y) -> h - basin(x, y),
        bcx=:open,
        bcy=:open,
        eta=(x, y) -> wave(x),
        u=(x, y) -> speed(x),
        dt=0.01 * res,
        n=0.016,
    )
end

function run(name="conicalisland", res=1)
    suffix = res == 1 ? "" : "_$res"
    outfile = to_path("out/$name$suffix.jld2")
    if isfile(outfile)
        @info "file exists: $(outfile)"
    else
        m = model(1)
        run!(m, outfile; seconds=15, frequency=25, output=(; eta=:eta))
        @info "file created: $(outfile)"
    end
end

#

function gauge_positions()
    d = 0.32 # depth (m)
    ad = 0.2 # a/d (case C)
    l = (2 * d / sqrt(0.75 * ad)) * acosh(sqrt(20))
    l2 = 12.96 - (3.6 + l / 2)
    gauge = [
        1 => (l2, 16.05, 0.320),
        2 => (l2, 14.55, 0.320),
        3 => (l2, 13.05, 0.320),
        4 => (l2, 11.55, 0.320),
        6 => (9.36, 13.80, 0.317),
        9 => (10.36, 13.80, 0.082),
        16 => (12.96, 11.22, 0.079),
        22 => (15.56, 13.80, 0.083),
    ]
    # shift to align centers (model - lab)
    sx, sy = (13, 15) .- (12.96, 13.80)
    Dict( g => (x=x + sx, y=y + sy, z=z) for (g, (x, y, z)) in gauge )
end

function plot_timeseries(name="conicalisland")
    gauge_ids = [2, 6, 9, 16, 22]
    gauge = gauge_positions()
    positions = [ (gauge[g].x, gauge[g].y) for g in gauge_ids ]

    # load data (laboratory)
    csv = CSV.File(to_path("in/conicalisland_timeseries_c.csv"); types=Float64)
    s = collect(csv.lookup[:s].column)
    lab = [ collect(csv.lookup[Symbol("g$g")].column) for g in gauge_ids ]

    # load data (model)
    filepath = to_path("out/$name.jld2")
    (; grid) = load(filepath, 0)
    sim = sample(filepath, gridindex.(Ref(grid), positions))

    # plot
    timeshift = 7.44 # s
    duration = 15.0 # s
    t1, tn = round.(Int, (timeshift, timeshift + duration) .* 25) .+ 1 # 25 hz
    r = t1:tn
    fig = Figure()
    header!(fig[1, 1], "Conical island, timeseries")
    limits = (s[r][1], s[r][end], -0.07, 0.12)
    rows = cld(length(gauge_ids), 2)
    for (i, (g, a, b)) in enumerate(zip(gauge_ids, lab, sim))
        ax = Axis(fig[mod1(i, rows) + 1, cld(i, rows)]; title="gauge $g", limits)
        lines!(ax, s[r], a[r]; color=:black, label="laboratory")
        lines!(ax, s[r], b[r .- (t1 - 1)]; color=:dodgerblue, label="model", linewidth=3)
        if i == rows + 1; axislegend(ax) end
    end
    fig
end

function plot_runup(name="conicalisland")
    frompolar(r, theta) = (cos(theta) * r, sin(theta) * r)
    function circle(center, r, n=60)
        section = 2 * pi / n
        [ frompolar(r, section * i) .+ center for i in 0:n ]
    end

    # load data (model)
    filepath = to_path("out/$name.jld2")
    (; grid, h) = load(filepath, 0)
    etamax = foldtime(filepath; init=zeros(grid.nx, grid.ny)) do eta, m
        max.(eta, m.eta)
    end

    # load data (laboratory)
    csv = CSV.File(to_path("in/conicalisland_runup_c.csv"); types=Float64)
    rad = collect(csv.lookup[:rad].column)
    runup = collect(csv.lookup[:runup].column) ./ 100 # vertical distance (from cm)
    # transform (from vertical to horizontal, from polar to cartesian)
    center = (13, 15)
    base = (3.6 - 1.1) * (0.305 / 0.625) # base distance between waterline and top
    waterline = 1.1 + base # radius of nominal waterline
    runup = waterline .- runup .* (base / 0.305) # distance from center
    runup = [ frompolar(r, t) .+ center for (r, t) in zip(runup, rad .- 0.5pi) ]

    # plot
    fig = Figure()
    header!(fig[1, 1], "Conical island, maximum run-up height")
    ax = Axis(fig[2, 1]; xlabel="x (m)", ylabel="y (m)",
        aspect=DataAspect(),
        limits=(8.5, 17.5, 11.5, 18.5),
        xticks=9:17,
        yticks=12:18
    )
    lines!(ax, circle(center, 1.1); color=:gray70) # cone top
    lines!(ax, circle(center, waterline); color=:gray70) # waterline
    scatter!(ax, runup; color=:black, markersize=6, label="laboratory")
    contour!(ax, grid.xc, grid.yc, (h .+ etamax);
        levels=[0.002], color=:dodgerblue, linewidth=3, label="model"
    )
    axislegend(ax)
    fig
end

function plot_results(name="conicalisland")
    filename = to_path("out/$name.jld2")
    plot_scene(Results(filename), "Conical island"; zscale=8.0, dz=0.05)
end

end # module
