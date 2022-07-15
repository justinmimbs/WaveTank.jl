
module SolitaryWave

using GLMakie
using Printf: @sprintf
using WaveTank: g, Grid, Model, run!, Results

import ..to_path
import ..plot_scene, ..plot_conservation!, ..header!

function solitarywave(a=0.1; h=1.0)
    c = sqrt(g * (h + a))
    k = sqrt(3a / 4h) / h
    l = 2pi / k
    eta(x) = a * sech(k * x)^2
    u(x) = let e = eta(x); c * e / (h + e) end
    (; a, h, c, k, l, eta, u)
end

function model(wave, nx)
    Model(;
        grid=Grid(
            (-1.5wave.l, 3.5wave.l), # x bounds
            (0.0, wave.l),           # y bounds
            (5nx, nx),               # grid resolution
        ),
        h=wave.h,                    # basin
        eta=(x, _) -> wave.eta(x),   # surface elevation
        u=(x, _) -> wave.u(x),       # particle velocity (x)
    )
end

function run(name="solitarywave"; ah=0.1, nx=100)
    outfile = to_path("out/$(name)_$(ah).jld2")
    if isfile(outfile)
        @info "file exists: $(outfile)"
    else
        wave = solitarywave(ah)
        period = wave.l / wave.c
        m = model(wave, nx)
        run!(m, outfile; seconds=2period, frequency=30 / period, output=(; eta=:eta))
        @info "file created: $(outfile)"
    end
end

function plot_comparison(name="solitarywave"; ah=0.1)
    filename = to_path("out/$(name)_$(ah).jld2")
    wave = solitarywave(ah)
    x = range(-1.5wave.l, 3.5wave.l; length=500)
    fig = Figure(; resolution=(1100, 600))
    header!(fig[1, 1], "Solitary wave")
    ax = Axis(fig[2, 1]; xlabel="x", ylabel="eta",)
    xlims!(ax, extrema(x))
    # analytical solution
    for i in 0:2
        lines!(ax, x, x -> wave.eta(x - i * wave.l); color=:red, label="solution", linewidth=3)
    end
    # model solution
    for (i, m) in enumerate(Results(filename))
        if mod1(i, 30) != 1; continue; end
        lines!(ax, m.grid.xc, m.eta[:, end ÷ 2]; color=:dodgerblue, label="model", linewidth=3)
    end
    axislegend(ax; merge=true)
    fig
end

function plot_conservation(name="solitarywave"; ah=0.1)
    filename = to_path("out/$(name)_$ah.jld2")
    fig = Figure()
    header!(fig[1, 1], "Solitary wave, mass conservation")
    plot_conservation!(fig[2, 1], Results(filename))
    fig
end

function plot_results(name="solitarywave"; ah=0.1)
    filename = to_path("out/$(name)_$ah.jld2")
    plot_scene(Results(filename), "Solitary wave"; zscale=20.0)
end

end # module
