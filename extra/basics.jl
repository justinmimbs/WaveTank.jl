using LinearAlgebra: norm
using WaveTank: g

# waves

function sech2wave(; a=1.0, k=1.0, x0=0.0)
    f(x) = a * sech(k * (x - x0))^2
end

function solitary_wave(a; h=1.0, x0=0.0)
    c = sqrt(g * (h + a))
    k = sqrt(3a / 4h) / h
    l = 2pi / k
    eta = sech2wave(; a, k, x0)
    u = particle_velocity(eta, h, a)
    (; a, h, c, k, l, eta, u)
end

# particle speed / phase speed = amplitude / depth
function particle_velocity(eta::Function, h, a)
    c = sqrt(g * (h + a))
    u(x) = let e = eta(x); c * e / (h + e) end
end


# basins

# 1D, waterline at x = 0
function planar_beach(depth, slope, slope_shore=slope)
    xs = [-(depth / slope), 0, 10depth / slope_shore]
    zs = [depth, 0, -10depth]
	piecewise_linear(xs, zs, depth)
end

function piecewise_linear(xs, ys, default=0.0)
    function f(x)
        if xs[1] <= x < xs[end]
            i = findfirst(>(x), xs) - 1
            x0, x1 = xs[i], xs[i + 1]
            y0, y1 = ys[i], ys[i + 1]
            alpha = (x - x0) / (x1 - x0)
            # lerp
            y0 + alpha * (y1 - y0)
        else
            default
        end
    end
end

function truncated_cone(x0, y0, rb, rt, h)
    l = rb - rt
    function f(x, y)
        alpha = (norm([x - x0, y - y0]) - rt) / l
        alpha <= 0.0 ? h : 1.0 <= alpha ? 0.0 : (1.0 - alpha) * h
    end
end
