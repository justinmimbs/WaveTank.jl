import Pkg
Pkg.activate(Base.current_project(@__DIR__))
cd(@__DIR__)

using GLMakie

using Neowave: Results
using Validation

function run()
    SolitaryWave.run()
    SimpleBeach.run(0.0185)
    SimpleBeach.run(0.3)
    ConicalIsland.run()
    MonaiValley.run()
end

function animate(obs::Observable, res::Results, seconds=Inf)
    for m in res
        if seconds <= m.dt * m.t; break; end
        obs[] = m
        sleep(0.016)
    end
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
# _, obs, res = MonaiValley.plot_results();
# animate(obs, res, 22.5)
