module WaveGeneration

using GLMakie
using Printf: @sprintf
using WaveTank: g, Grid, Model, run!, Results

import ..to_path
import ..sech2wave, ..particle_velocity
import ..plot_scene, ..plot_conservation!

function solitary_wave(ah; h=inv(ah), lx0=0.0)
    a = ah * h
    k = sqrt(3a / 4h) / h
    l = 8.3 / k
    c = sqrt(g * (h + a))
    eta = sech2wave(; a, x0=lx0 * l, k)
    u = particle_velocity(eta, a, h)
    (; a, h, k, l, c, eta, u)
end

function model(ah, gen, bcx)
    wave = solitary_wave(ah)
    if gen == :initial
        grid = Grid((-1.0wave.l, 3.0wave.l), (0.0, wave.l), (400, 100))
        eta(x, _) = wave.eta(x)
        u(x, _) = wave.u(x)
        m = Model(; grid, wave.h, bcx, eta, u)
        m, [], []
    else
        grid = Grid((1.0wave.l, 5.0wave.l), (0.0, wave.l), (400, 100))
        m = Model(; grid, wave.h, bcx)
        xt = range(1.0wave.l, -1.0wave.l; step=-(wave.c * m.dt))
        if gen == :wavemaker
            m, map(wave.u, xt), []
        else
            m, [], map(wave.eta, xt)
        end
    end
end

function run(name="solitarywave"; ah=0.1, gen=:initial, bcx=:wall)
    outfile = to_path("out/$(name)_$(ah)_$(gen)_$(bcx).jld2")
    if isfile(outfile)
        @info "file exists: $(outfile)"
    else
        wave = solitary_wave(ah)
        period = wave.l / wave.c
        seconds = period * (gen == :initial ? 2 : 4)
        m, wavemaker, waveinput = model(ah, gen, bcx)
        run!(m, outfile; seconds, frequency=inv(period), output=(; eta=:eta),
            wavemaker, waveinput,
        )
        @info "file created: $(outfile)"
    end
end

function runall(name="solitarywave")
    for gen in [:initial, :wavemaker, :waveinput], bcx in [:wall, :open]
        run(name; ah=0.1, gen, bcx)
    end
end

function plot_comparison!(gp, name="solitarywave"; ah=0.1, gen=:initial, bcx=:wall)
    results = Results(to_path("out/$(name)_$(ah)_$(gen)_$(bcx).jld2"))
    (; dt, grid) = results
    wave = solitary_wave(ah)

    # figure
    title = @sprintf("solitary wave: h = %d m, %s, %s", wave.h, gen, bcx)
    ax = Axis(gp; title, xlabel="x (m)", ylabel="elevation (m)")
    xlims!(ax, grid.xf[1], grid.xf[end])
    hidespines!(ax)

    for m in results
        # analytical solution
        lx = m.t * dt * wave.c
        lines!(ax, grid.xc, x -> wave.eta(x - lx); color=:gray)

        # model solution
        lines!(ax, grid.xc, m.eta[:, end ÷ 2]; color=:dodgerblue)
    end
    # TODO add legend; add time labels
    ax
end

function plot_comparison(name="solitarywave"; ah=0.1, gen=:initial, bcx=:wall)
    fig = Figure(resolution=(1000, 300))
    plot_comparison!(fig[1, 1], name; ah, gen, bcx)
    fig
end

function plot_comparisons(name="solitarywave";
    ah=[0.1],
    gen=[:initial, :wavemaker, :waveinput],
    bcx=[:wall, :open],
)
    args = [ (; ah, gen, bcx) for ah in ah for gen in gen for bcx in bcx ]
    fig = Figure(resolution=(1000, 240 * length(args)))
    for i in 1:length(args)
        plot_comparison!(fig[i, 1], name; args[i]...)
    end
    fig
end

function plot_conservation(name="solitarywave"; ah=0.1, gen=:initial, bcx=:wall)
    filename = to_path("out/$(name)_$(ah)_$(gen)_$(bcx).jld2")
    fig = Figure()
    plot_conservation!(fig[1, 1], Results(filename))
    fig
end

function plot_results(name="solitarywave"; ah=0.1, gen=:initial, bcx=:wall)
    filename = to_path("out/$(name)_$(ah)_$(gen)_$(bcx).jld2")
    plot_scene(Results(filename), "Solitary wave"; zscale=20.0)
end

end # module
