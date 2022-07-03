using Printf
using Neowave: Results

theme = let
    axis = (
        titlefont="TeX Gyre Heros",
        titlealign=:left,
    )
    Theme(
        resolution=(1100, 900),
        figure_padding=(20, 20, 20, 20),
        fontsize=18,
        Axis=(; axis...,
            leftspinecolor=:gray,
            rightspinecolor=:gray,
            bottomspinecolor=:gray,
            topspinecolor=:gray,
        ),
        Axis3=axis,
        Lines=(
            linewidth=2,
            color=:black,
        ),
        Colorbar=(
            spinewidth=0,
            ticklabelpad=5,
        ),
        Legend=(
            framecolor=:gray,
        ),
    )
end

function header!(gp::GridPosition, text)
    Label(gp, text; textsize=24, padding=(0, 0, 0, 0), halign=:left,
        tellwidth=false, tellheight=true
    )
end
function header!(fig::Figure, text)
    header!(fig[1, 1, TopLeft()], text)
end

function plot_conservation!(gp, iter)
    t_eta = [ m.t => sum(m.eta) for m in iter ]
    _, eta0 = t_eta[1]
    t = first.(t_eta)
    change = [ (eta / eta0 - 1.0) * 100 for (_, eta) in t_eta ]
    title = @sprintf("change in fluid: %.05f %%", change[end])
    ax = Axis(gp; title, xlabel="timestep", ylabel="change %")
    ylims!(ax, -0.5, 0.5)
    xlims!(ax, extrema(t))
    lines!(ax, t, change; color=:black)
end

function axis3(gp, m; zscale=5.0, zmax=nothing)
    h = m.h
    (; xf, yf) = m.grid
    zmin = -maximum(h)
    zmax = @something zmax max(-0.3zmin, -1.1minimum(h))
    Axis3(gp;
        aspect=(xf[end] - xf[1], yf[end] - yf[1], zscale * (zmax - zmin)),
        limits=(xf[1], xf[end], yf[1], yf[end], zmin, zmax),
        azimuth=1.33pi,
        elevation=0.25pi,
        xticks=WilkinsonTicks(6; k_min=3), xticklabelpad=2, xlabel="x (m)",
        yticks=WilkinsonTicks(5; k_min=3), yticklabelpad=2, ylabel="y (m)",
        zticks=WilkinsonTicks(3; k_min=3), zticklabelpad=2, zlabel="z (m)",
        protrusions=(60, 30, 30, 30)
    )
end

function plot_bathymetry!(ax, m; dz=0.0)
    b, xc, yc = -m.h, m.grid.xc, m.grid.yc
    zmin, zmax = extrema(b)
    dz = 0.0 < dz ? dz : zmax == zmin ? 1.0 : (zmax - zmin) / 21
    levels = vcat(-dz:-dz:zmin, dz:dz:zmax)
    contour3d!(ax, xc, yc, b; levels, color=:gray, overdraw=true)
    contour3d!(ax, xc, yc, b; levels=[0.0], color=:black, overdraw=true)
    sideprofiles = [
        (xc, fill(yc[1], length(xc)), b[:, 1]),
        (xc, fill(yc[end], length(xc)), b[:, end]),
        (fill(xc[1], length(yc)), yc, b[1, :]),
        (fill(xc[end], length(yc)), yc, b[end, :]),
    ]
    for (x, y, z) in sideprofiles
        lines!(ax, x, y, z; color=:gray, linewidth=1, overdraw=true)
    end
    ax
end

function plot_surface!(ax, m)
    (; grid, h, eta) = m
    etamax = 0.3 * maximum(h)
    depthmin = 0.005 * etamax
    dry = .<=(h .+ eta, depthmin)
    eta = copy(eta)
    eta[dry] .= NaN
    surface!(ax, grid.xc, grid.yc, eta;
        colormap=(Reverse(:balance), 0.9),
        colorrange=(-etamax, etamax),
        shading=false,
    )
end
function plot_surface!(ax, m::Observable)
    (; grid, h) = m[]
    etamax = 0.3 * maximum(h)
    depthmin = 0.005 * etamax
    surface = @lift let
        eta = copy(($m).eta)
        dry = .<=(h .+ eta, depthmin)
        eta[dry] .= NaN
        eta
    end
    surface!(ax, grid.xc, grid.yc, surface;
        colormap=(Reverse(:balance), 0.9),
        colorrange=(-etamax, etamax),
        shading=false,
    )
end

function plot_scene(m, title=""; zscale=5.0, zmax=nothing, dz=0.0)
    fig = Figure()
    header!(fig, title)
    ax = axis3(fig[1, 1], m; zscale, zmax)
    plot_bathymetry!(ax, m; dz)
    sp = plot_surface!(ax, m)
    Colorbar(fig[1, 2], sp)
    fig
end
function plot_scene(res::Results, title=""; zscale=5.0, zmax=nothing, dz=0.0)
    obs = Observable(first(res))
    fig = Figure()
    header!(fig, title)
    ax = axis3(fig[1, 1], obs[]; zscale, zmax)
    plot_bathymetry!(ax, obs[]; dz)
    sp = plot_surface!(ax, obs)
    Colorbar(fig[1, 2], sp)
    display(fig)
    fig, obs, res
end
