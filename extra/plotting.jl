function axis3(gp, m; zscale=5.0, zmax=nothing)
    h = m.h
    (; xf, yf) = m.grid
    zmin = -maximum(h)
    zmax = @something zmax max(-0.3zmin, -minimum(h))
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
    lines!(ax, xc, fill(yc[1], length(xc)), b[:, 1], color=:gray, overdraw=true)
    lines!(ax, xc, fill(yc[end], length(xc)), b[:, end], color=:gray, overdraw=true)
    lines!(ax, fill(xc[1], length(yc)), yc, b[1, :], color=:gray, overdraw=true)
    lines!(ax, fill(xc[end], length(yc)), yc, b[end, :], color=:gray, overdraw=true)
    contour3d!(ax, xc, yc, b; levels, color=:gray, overdraw=true)
    contour3d!(ax, xc, yc, b; levels=[0.0], color=:black, overdraw=true)
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

function plot_surface(m; zscale=5.0, zmax=nothing, dz=0.0)
    fig = Figure(resolution=(1000, 700))
    ax = axis3(fig[1, 1], m; zscale, zmax)
    plot_bathymetry!(ax, m; dz)
    sp = plot_surface!(ax, m)
    Colorbar(fig[1, 2], sp)
    fig
end
