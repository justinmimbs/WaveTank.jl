import Pkg
Pkg.activate(Base.current_project(@__DIR__))
cd(@__DIR__)

using GLMakie

using Neowave: Results
using Validation
using Validation: axis3, plot_bathymetry!, plot_surface!

function run()
    SolitaryWave.run()
    SimpleBeach.run(0.0185)
    SimpleBeach.run(0.3)
    ConicalIsland.run()
    MonaiValley.run()
end

# set_theme!(Validation.theme)

# SolitaryWave.plot_comparison()
# SolitaryWave.plot_conservation()

# SimpleBeach.plot_profiles()
# SimpleBeach.plot_conservation()

# ConicalIsland.plot_timeseries()
# ConicalIsland.plot_runup()

# MonaiValley.plot_setup()
# MonaiValley.plot_timeseries()

# function plot_obs(filename)
#     r = Results(filename)
#     m = Observable(first(r))
#     fig = Figure()
#     ax = axis3(fig[1, 1], m[]) # zscale, zmax
#     plot_bathymetry!(ax, m[]) # dz
#     sp = plot_surface!(ax, m)
#     Colorbar(fig[1, 2], sp)
#     display(fig)
#     fig, m, r
# end

# f1, m1, r1 = plot_obs("out/monaivalley.jld2");
# for m in r1
#     m1[] = m
#     sleep(0.016)
# end
