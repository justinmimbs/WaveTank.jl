import Pkg
Pkg.activate(Base.current_project(@__DIR__))

using BenchmarkTools
using ProfileView
using WaveTank

function model()
    h = 1.0
    a = 0.3
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

let
    m = step!(model())
    @info typeof(m.solver), hasproperty(m, :lu) ? typeof(m.lu) : nothing
end

function benchmark(samples=10)
    b1 = @benchmarkable step!(m, 200) setup=(m = model()) #teardown=(println(m.t))
    run(b1; evals=1, samples, seconds=300.0, gcsample=true)
end

function profile(steps=1)
    ProfileView.@profview step!(model(), steps)
end
