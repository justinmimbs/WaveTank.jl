module SimpleBeach

using CSV
using GLMakie
using Printf: @sprintf
using WaveTank

import ..to_path
import ..solitary_wave, ..planar_beach
import ..plot_scene, ..plot_conservation!, ..header!

function model(ah=0.3) # ah = amplitude / depth ratio
    h = 1.0
    a = ah * h
    k = sqrt(3 * a / (4 * h^3))
    l = (2.0 / k) * acosh(sqrt(0.05^-1))
    slope = 1 / 19.85
    bathymetry = planar_beach(h, slope)
    wave = solitary_wave(a; h, x0=-(0.5l + h / slope))
    grid = let
        dx = h / 8
        xmin = -(1.5l + h / slope)
        nx = ceil(Int, (1.5l + 2h / slope) / dx)
        xmax = xmin + nx * dx
        Grid((xmin, xmax), (0, 80dx), (nx, 80))
    end
    Model(;
        grid,
        h=(x, y) -> bathymetry(x),
        bcx=:wall,
        bcy=:wall,
        eta=(x, y) -> wave.eta(x),
        u=(x, y) -> wave.u(x),
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
    results = Results(to_path("out/$(name)_$ah.jld2"))
    ts = ah == 0.3 ? (15:5:30) : (30:10:70)
    limits = ah == 0.3 ? (-20, 8, -0.05, 0.45) : (-20, 3.5, -0.03, 0.08)

    fig = Figure()
    header!(fig[1, 1], "Simple beach, profiles")
    for (i, t) in enumerate(ts)
        (; grid, h, eta) = results[32 * t - 20] # s = 0.32 * t - 0.2
        j = length(grid.yc) ÷ 2
        ax = Axis(fig[1 + i, 1]; title="t = $t", limits)
        scatter!(ax, -csv["x_t$t"], csv["eta_t$t"]; color=:black, markersize=5, label="laboratory")
        lines!(ax, grid.xc, eta[:, j]; color=:dodgerblue, linewidth=3, label="model")
        lines!(ax, grid.xc, -h[:, j]; color=:gray, linewidth=2)
        if i == 1
            axislegend(ax)
        end
    end
    fig
end

function plot_conservation(name="simplebeach"; ah=0.3)
    filename = to_path("out/$(name)_$ah.jld2")
    fig = Figure()
    header!(fig[1, 1], "Simple beach, mass conservation")
    plot_conservation!(fig[2, 1], Results(filename))
    fig
end

function plot_results(name="simplebeach"; ah=0.3)
    filename = to_path("out/$(name)_$ah.jld2")
    plot_scene(Results(filename), "Simple beach"; zscale=ah == 0.3 ? 5.0 : 20.0)
end

end # module
