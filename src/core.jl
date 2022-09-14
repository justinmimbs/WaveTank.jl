using OffsetArrays: OffsetArrays, OffsetArray

const g = 9.80665 # m s^-2

grid_doc = """
    Grid(xbounds, ybounds, resolution)

Define a 2D computational grid by its domain bounds and resolution in each dimension.

```jldoctest
julia> Grid((-15.0, 10.0), (-5.0, 5.0), (500, 200))
Grid 500 × 200 (25 m × 10 m)
    Δx = 0.05 m, bounds x = [-15, 10]
    Δy = 0.05 m, bounds y = [-5, 5]
```
"""
struct Grid
    nx::Int # number of cells
    ny::Int # number of cells
    dx::Float64 # cell length x
    dy::Float64 # cell length y
    xc::Vector{Float64} # center values (length = nx)
    xf::Vector{Float64} # face values (length = nx + 1)
    yc::Vector{Float64} # center values (length = ny)
    yf::Vector{Float64} # face values (length = ny + 1)
    @doc grid_doc
    function Grid(x, y, (nx, ny))
        lx, ly = x[end] - x[1], y[end] - y[1]
        dx, dy = lx / nx, ly / ny
        xc = range(x[1], x[end] - dx, length=nx) .+ 0.5dx
        xf = range(x[1], x[end], length=nx + 1)
        yc = range(y[1], y[end] - dy, length=ny) .+ 0.5dy
        yf = range(y[1], y[end], length=ny + 1)
        new(
            nx, ny, dx, dy,
            collect(xc), collect(xf), collect(yc), collect(yf),
        )
    end
end

function Base.show(io::IO, grid::Grid)
    (; nx, ny, xf, yf) = grid
    @printf(io, "Grid %d × %d (%g m × %g m)", nx, ny, xf[end] - xf[1], yf[end] - yf[1])
end
function Base.show(io::IO, ::MIME"text/plain", grid::Grid)
    (; nx, ny, xf, yf, dx, dy) = grid
    xl, xu = xf[1], xf[end]
    yl, yu = yf[1], yf[end]
    @printf(io, "Grid %d × %d (%g m × %g m)", nx, ny, xu - xl, yu - yl)
    @printf(io, "\n    Δx = %g m, bounds x = [%g, %g]", dx, xl, xu)
    @printf(io, "\n    Δy = %g m, bounds y = [%g, %g]", dy, yl, yu)
end

# array utilities

function interior(v::Vector)
    @view v[2:end - 1]
end
function interior(a::Matrix)
    @view a[2:end - 1, 2:end - 1]
end
function interiorx(a)
    @view a[2:end - 1, :]
end
function interiory(a)
    @view a[:, 2:end - 1]
end

function diffx(a)
    @views a[2:end, :] .- a[1:end - 1, :]
end
function diffy(a)
    @views a[:, 2:end] .- a[:, 1:end - 1]
end

function averagex(a)
    @views 0.5 .* (a[1:end - 1, :] .+ a[2:end, :])
end
function averagey(a)
    @views 0.5 .* (a[:, 1:end - 1] .+ a[:, 2:end])
end
function average(a)
    @inbounds [ 0.25 * (a[i, j] + a[i, j + 1] + a[i + 1, j] + a[i + 1, j + 1])
        for i in 1:size(a, 1) - 1, j in 1:size(a, 2) - 1
    ]
end

function min!(a, v)
    @inbounds for i in 1:length(a); a[i] = min(a[i], v); end
    a
end
function max!(a, v)
    @inbounds for i in 1:length(a); a[i] = max(a[i], v); end
    a
end

# ghost cells by linear extrapolation
# w = 0: equivalant to zero Neumann condition (copies neighboring values)
# w = 1: equivalent to Neumann condition using change from neighbors (extrapolates values)
function halo(a::Matrix, w=0.0)
    b = zeros(size(a) .+ (2, 2))
    interior(b) .= a
    haloy!(halox!(b, w), w)
end
function halox(a, w=0.0)
    b = zeros(size(a) .+ (2, 0))
    interiorx(b) .= a
    halox!(b, w)
end
function haloy(a, w=0.0)
    b = zeros(size(a) .+ (0, 2))
    interiory(b) .= a
    haloy!(b, w)
end
function halox!(b, w=0.0)
    if w == 0.0
        @views b[1, :] .= b[2, :]
        @views b[end, :] .= b[end - 1, :]
    else
        @views b[1, :] .= b[2, :] .- w .* (b[3, :] .- b[2, :]) # x-
        @views b[end, :] .= b[end - 1, :] .- w .* (b[end - 2, :] .- b[end - 1, :]) # x+
    end
    b
end
function haloy!(b, w=0.0)
    if w == 0.0
        @views b[:, 1] .= b[:, 2]
        @views b[:, end] .= b[:, end - 1]
    else
        @views b[:, 1] .= b[:, 2] .- w .* (b[:, 3] .- b[:, 2]) # y-
        @views b[:, end] .= b[:, end - 1] .- w .* (b[:, end - 2] .- b[:, end - 1]) # y+
    end
    b
end


# depthsmooth

function extend(a, k::Int)
    nx, ny = size(a)
    b = OffsetArray(zeros(nx + 2k, ny + 2k), -k, -k)
    b[1:nx, 1:ny] .= a
    for i in  1 .- (1:k); b[i, :] .= b[ 1, :]; end
    for i in nx .+ (1:k); b[i, :] .= b[nx, :]; end
    for j in  1 .- (1:k); b[:, j] .= b[:,  1]; end
    for j in ny .+ (1:k); b[:, j] .= b[:, ny]; end
    b
end
# Create a normalized kernel, scaled so that r spans 3 sigmas, and the kernel size is an odd number.
function gaussiankernel(r::Float64)
    s = (3.0 / r)^2 / 2.0
    n = max(1, round(Int, r) - 1)
    k = [ exp(-(i^2 + j^2) * s) for i in -n:n, j in -n:n ]
    k .*= inv(sum(k))
    OffsetArrays.centered(k)
end
function depthsmooth(h, alpha) # expect alpha scaled by dx
    hr = h * alpha * 1.5 # search radius at each location, 1.5 = sqrt(9/4)
    hx = extend(h, ceil(Int, maximum(hr))) # h extended
    Base.map(CartesianIndices(h)) do i
        r = hr[i]
        if r < 0.5
            hx[i]
        else
            k = gaussiankernel(r)
            ki = CartesianIndices(k)
            sum((@view hx[ki .+ i]) .* k)
        end
    end
end


# advection

# split a velocity component into two fields, for positive and negative flows
function upwind(a)
    # p, n: faces (14)
    max.(a, 0.0), min.(a, 0.0)
end

# split u (momentum-conserving)
# ∘|⋅|⋅|∘ where ∘ = (eta); | = (u); ⋅ = (h, d)
function upwindu(u, eta, h, d)
    nx, ny = size(h)
    # centers
    umean = averagex(u)
    # left flux: rightmost interior faces (22a)
    etap = @inbounds ( 0.0 < u[i - 1, j] ? (eta[i - 1, j] + eta[i, j]) / 2 : eta[i, j]
        for i in 2:nx + 1, j in 1:ny
    )
    flup = similar(u)
    @views flup[2:end, :] .= umean .* (h .+ etap)
    # right flux: leftmost interior faces (22b)
    etan = @inbounds ( u[i + 1, j] < 0.0 ? (eta[i + 2, j] + eta[i + 1, j]) / 2 : eta[i + 1, j]
        for i in 1:nx, j in 1:ny
    )
    flun = similar(u)
    @views flun[1:end - 1, :] .= umean .* (h .+ etan)
    # uhat: interior faces (21)
    dmean = averagex(d)
    @views flup[2:end - 1, :] ./= dmean
    @views flun[2:end - 1, :] ./= dmean
    # up, un: faces (20)
    halox!(max!(flup, 0.0)), halox!(min!(flun, 0.0))
end

# split v (momentum-conserving)
# ∘|⋅|⋅|∘ transposed, where ∘ = (eta); | = (v); ⋅ = (h, d)
function upwindv(v, eta, h, d)
    nx, ny = size(h)
    # centers
    vmean = averagey(v)
    # pos flux: endmost interior faces (25a)
    etap = @inbounds ( 0.0 < v[i, j - 1] ? (eta[i, j - 1] + eta[i, j]) / 2 : eta[i, j]
        for i in 1:nx, j in 2:ny + 1
    )
    flvp = similar(v)
    @views flvp[:, 2:end] .= vmean .* (h .+ etap)
    # neg flux: firstmost interior faces (25b)
    etan = @inbounds ( v[i, j + 1] < 0.0 ? (eta[i, j + 2] + eta[i, j + 1]) / 2 : eta[i, j + 1]
        for i in 1:nx, j in 1:ny
    )
    flvn = similar(v)
    @views flvn[:, 1:end - 1] .= vmean .* (h .+ etan)
    # vhat: interior faces (24)
    dmean = averagey(d)
    @views flvp[:, 2:end - 1] ./= dmean
    @views flvn[:, 2:end - 1] ./= dmean
    # vp, vn: faces (23)
    haloy!(max!(flvp, 0.0)), haloy!(min!(flvn, 0.0))
end
