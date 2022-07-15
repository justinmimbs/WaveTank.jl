import Pkg
Pkg.activate(Base.current_project(@__DIR__))
cd(@__DIR__)

using GLMakie

using Validation
using Validation: animate

function run()
    SolitaryWave.run()
    SimpleBeach.run(0.0185)
    SimpleBeach.run(0.3)
    ConicalIsland.run()
    MonaiValley.run()
end

set_theme!(Validation.theme)

#=
SolitaryWave.plot_comparison()
SolitaryWave.plot_conservation()

SimpleBeach.plot_profiles(ah=0.0185)
SimpleBeach.plot_conservation(ah=0.0185)
animate(SimpleBeach.plot_results(ah=0.3))

SimpleBeach.plot_profiles(ah=0.3)
SimpleBeach.plot_conservation(ah=0.3)
animate(SimpleBeach.plot_results(ah=0.0185))

ConicalIsland.plot_setup()
ConicalIsland.plot_timeseries()
ConicalIsland.plot_runup()
animate(ConicalIsland.plot_results())

MonaiValley.plot_setup()
MonaiValley.plot_timeseries()
animate(MonaiValley.plot_results(), 22.5)
=#
