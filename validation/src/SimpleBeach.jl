module SimpleBeach

using CSV
using GLMakie
using Printf: @sprintf
using Neowave: g, Grid, Model, run!, Results, load

import ..to_path
import ..sech2wave, ..particle_velocity, ..piecewiselinear
import ..plot_scene, ..plot_conservation!

function model(ah=0.3) # ah = amplitude / depth ratio
    h = 1.0
    a = ah * h
    k = sqrt(3 * a / (4 * h^3))
    l = (2.0 / k) * acosh(sqrt(0.05^-1))
    beta = 1 / 19.85 # beach slope
    # beach (waterline at x = 0)
    bmin = -h / beta
    bh = 2 * h # beach height
    bmax = bmin + (bh / beta)
    # x domain
    xmin = bmin - 1.5 * l
    dx = h / 8
    nx = ceil(Int, (bmax - xmin) / dx)
    xmax = xmin + nx * dx
    #
    wave = sech2wave(a=a, x0=bmin - 0.5 * l, k=k)
    speed = particle_velocity(wave, a, h)
    bathymetry = piecewiselinear([bmin, bmax, bmax + l], h .- [0, bh, bh], h)
    Model(
        grid=Grid((xmin, xmax), (0, 80dx), (nx, 80)),
        h=(x, y) -> bathymetry(x),
        bcx=:wall,
        bcy=:wall,
        eta=(x, y) -> wave(x),
        u=(x, y) -> speed(x),
        n=0.01,
        dt=0.01,
    )
end

function run(name="simplebeach"; ah=0.3)
    outfile = to_path("out/$(name)_$ah.jld2")
    if isfile(outfile)
        @info "file exists: $(outfile)"
    else
        m = model(ah)
        run!(m, outfile; seconds=30, frequency=5, output=(; eta=:eta))
        @info "file created: $(outfile)"
    end
end

function plot_profiles(name="simplebeach"; ah=0.3)
    csv = CSV.File(to_path("in/simplebeach_profiles_$ah.csv"), types=Float64)
    ts = ah == 0.3 ? (15:5:30) : (30:10:70)
    limits = ah == 0.3 ? (-8, 20, -0.05, 0.45) : (-2, 20, -0.03, 0.08)

    fig = Figure(; resolution=(800, 240 * length(ts)))
    for (i, t) in enumerate(ts)
        (; grid, h, eta) = load(to_path("out/$(name)_$ah.jld2"), 32 * t - 20) # s = 0.32 * t - 0.2
        j = length(grid.yc) ÷ 2
        ax = Axis(fig[i, 1]; title="t = $t", limits)
        scatter!(ax, csv["x_t$t"], csv["eta_t$t"]; color=:black, markersize=3)
        lines!(ax, -reverse(grid.xc), reverse(eta[:, j]); color=:dodgerblue, linewidth=2)
        lines!(ax, -reverse(grid.xc), -reverse(h[:, j]); color=:gray, linewidth=2)
    end
    fig
end

function plot_conservation(name="simplebeach"; ah=0.3)
    filename = to_path("out/$(name)_$ah.jld2")
    fig = Figure()
    plot_conservation!(fig[1, 1], Results(filename))
    fig
end

function plot_results(name="simplebeach"; ah=0.3)
    filename = to_path("out/$(name)_$ah.jld2")
    plot_scene(Results(filename), "Simple beach"; zscale=5.0)
end

end # module
