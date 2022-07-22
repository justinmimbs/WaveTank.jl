using LinearAlgebra: norm
using WaveTank: g

function sech2wave(; a=1.0, x0=0.0, k=1)
    # a : amplitude, k : wavenumber
    f(x) = a * sech(k * (x - x0))^2
end

# particle speed / phase speed = amplitude / depth
function particle_velocity(eta::Function, h, a)
    c = sqrt(g * (h + a))
    u(x) = let e = eta(x); c * e / (h + e) end
end

function piecewiselinear(xs, ys, default=0.0)
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

function truncatedcone(x0, y0, rb, rt, h)
    l = rb - rt
    function f(x, y)
        alpha = (norm([x - x0, y - y0]) - rt) / l
        alpha <= 0.0 ? h : 1.0 <= alpha ? 0.0 : (1.0 - alpha) * h
    end
end
