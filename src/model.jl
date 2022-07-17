using LinearAlgebra: LinearAlgebra, Adjoint, norm, mul!, ldiv!
using SparseArrays: SparseMatrixCSC, spdiagm, nonzeros
using ILUZero: ILU0Precon, ilu0, ilu0!
using ConjugateGradients: BiCGStabData, bicgstab!

mutable struct Model
    # static
    grid::Grid
    bcx::Symbol # boundary conditions x-normals
    bcy::Symbol # boundary conditions y-normals
    n::Float64 # Manning roughness coefficient
    dt::Float64
    h::Matrix{Float64} # centers
    hi::Matrix{Float64} # intersections
    hfx::Matrix{Float64} # faces x
    hfy::Matrix{Float64} # faces y
    h_x::Matrix{Float64} # faces x
    h_y::Matrix{Float64} # faces y
    # time-stepped
    t::Int
    eta::Matrix{Float64} # centers
    u::Matrix{Float64} # faces x
    v::Matrix{Float64} # faces y
    ws::Matrix{Float64} # centers
    wb::Matrix{Float64} # centers
    # workspace
    coeff::SparseMatrixCSC{Float64, Int}
    lu::ILU0Precon{Float64, Int}
    solver::BiCGStabData{Float64}
    q::Matrix{Float64} # centers
    #
    function Model(; grid, h=1.0, u=nothing, v=nothing, eta=nothing,
        bcx=:wall, bcy=:wall, n=0.0, dt=nothing, courant=0.2,
        hsmooth=0.0,
    )
        @assert bcx in [:wall, :open]
        @assert bcy in [:wall, :open]
        grid = convert(Grid, grid)
        nx, ny = grid.nx, grid.ny
        xf, yf, xc, yc = grid.xf, grid.yf, grid.xc, grid.yc
        # essential
        h = isa(h, Number) ? fill(h, nx, ny) : [ h(x, y) for x in xc, y in yc ]
        if 0.0 < hsmooth
            h = depthsmooth(h, hsmooth / grid.dx)
        end
        u = u === nothing ? zeros(nx + 1, ny) : [ u(x, y) for x in xf, y in yc ]
        v = v === nothing ? zeros(nx, ny + 1) : [ v(x, y) for x in xc, y in yf ]
        eta = eta === nothing ? zeros(nx, ny) : [ eta(x, y) for x in xc, y in yc ]
        eta .= max.(eta, -h) # ensure 0 <= depth
        # derived
        hi = average(h)
        hhalo = halo(h, 1.0)
        hfx = averagex(interiory(hhalo))
        hfy = averagey(interiorx(hhalo))
        h_x = diffx(interiory(hhalo))
        h_y = diffy(interiorx(hhalo))
        ws = -(eta .+ h) .* (diffx(u) ./ grid.dx .+ diffy(v) ./ grid.dy)
        wb = -(averagex(u) .* diffx(hfx) ./ grid.dx) .+
            -(averagey(v) .* diffy(hfy) ./ grid.dy)
        # workspace
        coeff = let
            # pre-allocate diagonals
            diag1 = zeros((nx * ny) - 1)
            diagnx = zeros(nx * (ny - 1))
            spdiagm(
                -nx => diagnx,
                -1 => diag1,
                0 => zeros(Float64, nx * ny),
                1 => diag1,
                nx => diagnx,
            )
        end
        lu = ilu0(coeff)
        solver = BiCGStabData(nx * ny, Float64)
        q = zeros(nx, ny)
        # dt (if not given, then set by Courant number)
        if !isa(dt, Number)
            cg = sqrt(g * maximum(h + eta)) # group speed
            dt = courant * min(grid.dx, grid.dy) / cg
            dt = inv(ceil(inv(dt))) # ensure hz is an integer (i.e. dt = 1/{integer})
        end
        #
        new(grid, bcx, bcy, n, dt, h, hi, hfx, hfy, h_x, h_y,
            0, eta, u, v, ws, wb,
            coeff, lu, solver, q
        )
    end
end

function Base.show(io::IO, m::Model)
    (; t, dt) = m
    print(io, "WaveTank.Model @ t = $t ($(t * dt) s)")
end
function Base.show(io::IO, ::MIME"text/plain", m::Model)
    (; t, dt, grid, bcx, bcy) = m
    println(io, "WaveTank.Model @ t = $t ($(t * dt) s)")
    println(io, "   grid = $grid")
    println(io, "   boundaries = (x = $bcx, y = $bcy)")
    print(  io, "   Δt = $dt, Δx = $(grid.dx), Δy = $(grid.dy)")
end

function step!(m::Model)::Model
    dt, dx, dy = m.dt, m.grid.dx, m.grid.dy
    nx, ny = m.grid.nx, m.grid.ny

    # centers
    etahalo = halo(m.eta)
    etahalox = interiory(etahalo)
    etahaloy = interiorx(etahalo)

    # 1. wet/dry evaluation
    # update! m.eta
    d, dry, dryfx, dryfy = @views let
        # d, wet: centers
        d = m.eta + m.h
        wet = [ 100 * eps() < d for d in d ]

        # advance waterline (into dry cells)
        etanext = copy(m.eta)
        inflow::Vector{Float64} = []
        for j in 1:ny, i in 1:nx
            # skip wet cells
            if wet[i, j]
                continue
            end
            # if water is inflowing, store change
            # from x-
            if 1 < i && wet[i - 1, j] #&& 0 < m.u[i, j]
                eta = 0.0 - m.u[i, j] * d[i - 1, j] / dx
                push!(inflow, eta)
            end
            # from x+
            if i < nx && wet[i + 1, j] #&& m.u[i + 1, j] < 0
                eta = m.u[i + 1, j] * d[i + 1, j] / dx
                push!(inflow, eta)
            end
            # from y-
            if 1 < j && wet[i, j - 1] #&& 0 < m.v[i, j]
                eta = 0 - m.v[i, j] * d[i, j - 1] / dy
                push!(inflow, eta)
            end
            # from y+
            if j < ny && wet[i, j + 1] #&& m.v[i, j + 1] < 0
                eta = m.v[i, j + 1] * d[i, j + 1] / dy
                push!(inflow, eta)
            end
            # if any inflow: update eta with change
            if !isempty(inflow)
                etanext[i, j] = m.eta[i, j] - dt * sum(inflow)
                empty!(inflow)
            end
        end
        m.eta = etanext

        d = m.eta + m.h
        wet = [ 100 * eps() < d for d in d ]
        dry = .!wet

        # wetfx, dryfx: x faces
        wetfx = fill(true, nx + 1, ny)
        wetfx[1:end - 1, :] .= wet
        wetfx[2:end, :] .&= wet
        dryfx = .!wetfx

        # wetfx, dryfx: x faces
        wetfy = fill(true, nx, ny + 1)
        wetfy[:, 1:end - 1] .= wet
        wetfy[:, 2:end] .&= wet
        dryfy = .!wetfy

        # ensure non-zero depth
        d[dry] .= 100 * eps()
        m.eta[dry] .= -m.h[dry] # new!
        # domain bounds
        m.u[dryfx] .= 0.0
        m.v[dryfy] .= 0.0

        d, dry, dryfx, dryfy
    end

    # 2. velocity advection (approximate)
    # u, v

    di = average(d)

    # u: x faces (18)
    u = let
        # eta_x: x faces, interior
        eta_x = diffx(m.eta) ./= dx
        # up, un: x faces, interior
        up, un = map(interiorx, upwindu(m.u, etahalox, m.h, d))
        # u_x: centers
        u_x = diffx(m.u) ./= dx
        # uu_x: x faces, interior (u * upwind u_x)
        # uu_x = @views (up .* u_x[1:end - 1, :] .+ un .* u_x[2:end, :])
        uu_x = up .*= @view u_x[1:end - 1, :]
        uu_x .+= un .* @view u_x[2:end, :]
        # vx, vxp, vxn: x faces, interior
        vx = average(m.v)
        vxp, vxn = upwindv(vx, average(etahaloy), m.hi, di)
        # u_y: corners, x interior
        u_y = diffy(haloy(interiorx(m.u))) ./= dy
        # vu_y: x faces, interior (v * upwind u_y)
        # vu_y = @views (vxp .* u_y[:, 1:end - 1] .+ vxn .* u_y[:, 2:end])
        vu_y = vxp .*= @view u_y[:, 1:end - 1]
        vu_y .+= vxn .* @view u_y[:, 2:end]
        # friction term, x interior
        fr = let
            s = m.n^2 * g
            ux = interiorx(m.u)
            a = averagex(d)
            @inbounds for i in 1:length(a)
                a[i] = s * ux[i] * norm((ux[i], vx[i])) / a[i]^(4.0/3.0)
            end
            a
        end
        # u: x faces (interior update)
        u = copy(m.u)
        interiorx(u) .-= m.dt .* (g .* eta_x .+ uu_x .+ vu_y .+ fr)
        # domain bounds
        u[dryfx] .= 0.0
        u
    end

    # v: y faces (19)
    v = let
        # eta_y: y faces, interior
        eta_y = diffy(m.eta) ./= dy
        # vp, vn: y faces, interior
        vp, vn = map(interiory, upwindv(m.v, etahaloy, m.h, d))
        # v_y: centers
        v_y = diffy(m.v) ./= dy
        # vv_y: y faces, interior (v * upwind v_y)
        # vv_y = @views (vp .* v_y[:, 1:end - 1] .+ vn .* v_y[:, 2:end])
        vv_y = vp .*= @view v_y[:, 1:end - 1]
        vv_y .+= vn .* @view v_y[:, 2:end]
        # uy, uyp, uyn: y faces, interior
        uy = average(m.u)
        uyp, uyn = upwindu(uy, average(etahalox), m.hi, di)
        # v_x: corners, y interior
        v_x = diffx(halox(interiory(m.v))) ./= dx
        # uv_x: y faces, interior (u * upwind v_x)
        # uv_x = @views (uyp .* v_x[1:end - 1, :] .+ uyn .* v_x[2:end, :])
        uv_x = uyp .*= @view v_x[1:end - 1, :]
        uv_x .+= uyn .* @view v_x[2:end, :]
        # friction term, y interior
        fr = let
            s = m.n^2 * g
            vy = interiory(m.v)
            a = averagey(d)
            @inbounds for i in 1:length(a)
                a[i] = s * vy[i] * norm((uy[i], vy[i])) / a[i]^(4.0/3.0)
            end
            a
        end
        # v: y faces (interior update)
        v = copy(m.v)
        interiory(v) .-= m.dt .* (g .* eta_y .+ uv_x .+ vv_y .+ fr)
        # domain bounds
        v[dryfy] .= 0.0
        v
    end

    # 3. pressure correction
    # update! u, m.wb, m.ws, m.q, m.coeff
    let
        # a, b: x faces, y faces (31)
        a, b = let
            etamh = m.eta - m.h
            a = @views (diffx(etamh) ./ (d[1:end - 1, :] .+ d[2:end, :]))
            b = @views (diffy(etamh) ./ (d[:, 1:end - 1] .+ d[:, 2:end]))
            halox(a), haloy(b)
        end
        # a[dryfx] .= 0.0
        # b[dryfy] .= 0.0

        # coeff: centers (38)
        let
            sx = m.dt / (2 * dx^2)
            sy = m.dt / (2 * dy^2)
            dt2 = 2 * m.dt
            coeff = m.coeff
            fill!(nonzeros(coeff), 0.0)
            # build transposed (i.e. set coeff[j, i] instead of [i, j])
            @inbounds for j in 1:ny, i in 1:nx
                ij = (j - 1) * nx + i
                # center
                center = dt2 / d[i, j]^2
                # x- (west)
                if 1 < i && !dry[i - 1, j]
                    coeff[ij - 1, ij] = sx * (-1 + a[i, j])
                    center           += sx * ( 1 + a[i, j])
                end
                # x+ (east)
                if i < nx && !dry[i + 1, j]
                    coeff[ij + 1, ij] = sx * (-1 - a[i + 1, j])
                    center           += sx * ( 1 - a[i + 1, j])
                end
                # y- (south)
                if 1 < j && !dry[i, j - 1]
                    coeff[ij - nx, ij] = sy * (-1 + b[i, j])
                    center            += sy * ( 1 + b[i, j])
                end
                # y+ (north)
                if j < ny && !dry[i, j + 1]
                    coeff[ij + nx, ij] = sy * (-1 - b[i, j + 1])
                    center            += sy * ( 1 - b[i, j + 1])
                end
                coeff[ij, ij] = center
            end
        end

        # wb: centers (33)
        wb = let
            umean = averagex(m.u) # centers
            vmean = averagey(m.v) # centers
            @inbounds [
                -umean[i, j] * m.h_x[0 < umean[i, j] ? i : i + 1, j] +
                -vmean[i, j] * m.h_y[i, 0 < vmean[i, j] ? j : j + 1]
                for i in 1:nx, j in 1:ny
            ]
        end

        # div: centers (39)
        div = zeros(nx, ny)
        div .-= diffx(u) ./ dx .+ diffy(v) ./ dy .+ (m.ws .+ m.wb .- 2 .* wb) ./ d
        # div[dry] .= 0.0

        # q: centers
        let
            ilu0!(m.lu, m.coeff)
            precon(out, x) = ldiv!(out, m.lu', x)
            coeff(out, x) = mul!(out, m.coeff', x)
            bicgstab!(coeff, vec(div), vec(m.q); data=m.solver, precon, tol=1e-08)
        end
        m.q[dry] .= 0.0

        # ws: centers (32)
        m.ws .+= m.wb .- wb .+ (2 * m.dt) .* m.q ./ d
        m.wb = wb

        # u, update: faces, interior (29)
        interiorx(u) .-= (dt / dx) .* (interiorx(a) .* averagex(m.q) .+ diffx(m.q) ./ 2)

        # v, update: faces, interior (30)
        interiory(v) .-= (dt / dy) .* (interiory(b) .* averagey(m.q) .+ diffy(m.q) ./ 2)
    end

    # speed check: excessive velocities are probably spurious, so reduce them
    @inbounds let
        cmax = 0.5 * min(dx, dy) / dt
        cmin = sqrt(eps())
        for i in 1:length(u)
            u[i] = u[i] < -cmax ? -cmin : cmax < u[i] ? cmin : u[i]
        end
        for i in 1:length(v)
            v[i] = v[i] < -cmax ? -cmin : cmax < v[i] ? cmin : v[i]
        end
    end

    # 4. boundary conditions part 1
    # half-step (0.5dt)
    # d @ (t-1, x-1); u @ (t-1, x-1); no soliton artifacts, but unstable edges
    if m.bcx == :open
        # u (x normal)
        d = max.(m.eta .+ m.h, 0.0)
        c1x, cnx = @views (sqrt.(g .* d[2, :]), sqrt.(g .* d[end - 1, :]))
        @views u[1, :] .-= 0.5dt .* (u[1, :] .- c1x) .* (m.u[3, :] .- m.u[2, :]) ./ dx
        @views u[end, :] .-= 0.5dt .* (u[end, :] .+ cnx) .* (m.u[end - 1, :] .- m.u[end - 2, :]) ./ dx
    end
    if m.bcy == :open
        # v (y normal)
        d = max.(m.eta .+ m.h, 0.0)
        c1y, cny = @views (sqrt.(g .* d[:, 2]), sqrt.(g .* d[:, end - 1]))
        @views v[:, 1] .-= 0.5dt .* (v[:, 1] .- c1y) .* (m.v[:, 3] .- m.v[:, 2]) ./ dy
        @views v[:, end] .-= 0.5dt .* (v[:, end] .+ cny) .* (m.v[:, end - 1] .- m.v[:, end - 2]) ./ dy
    end

    # 5. free-surface extraction
    # update! m.eta
    let
        # flx: faces (13a)
        flx = let
            up, un = upwind(u)
            ueta = @views (up .* etahalox[1:end - 1, :] .+ un .* etahalox[2:end, :])
            ueta .+= (u .* m.hfx)
        end
        # fly: faces (13b)
        fly = let
            vp, vn = upwind(v)
            veta = @views (vp .* etahaloy[:, 1:end - 1] .+ vn .* etahaloy[:, 2:end])
            veta .+= (v .* m.hfy)
        end
        # eta: centers (12)
        m.eta .-= m.dt .* (diffx(flx) ./ dx .+ diffy(fly) ./ dy)
    end

    # 6. boundary conditions part 2
    # half-step (0.5dt)
    # d @ (t, x); u @ (t-1, x); stable edges, but soliton artifacts
    if m.bcx == :open
        # u (x normal)
        d = max.(m.eta .+ m.h, 0.0)
        c1x, cnx = @views (sqrt.(g .* d[1, :]), sqrt.(g .* d[end, :]))
        @views u[1, :] .-= 0.5dt .* (m.u[1, :] .- c1x) .* (m.u[2, :] .- m.u[1, :]) ./ dx
        @views u[end, :] .-= 0.5dt .* (m.u[end, :] .+ cnx) .* (m.u[end, :] .- m.u[end - 1, :]) ./ dx
    end
    if m.bcy == :open
        # v (y normal)
        d = max.(m.eta .+ m.h, 0.0)
        c1y, cny = @views (sqrt.(g .* d[:, 1]), sqrt.(g .* d[:, end]))
        @views v[:, 1] .-= 0.5dt .* (m.v[:, 1] .- c1y) .* (m.v[:, 2] .- m.v[:, 1]) ./ dy
        @views v[:, end] .-= 0.5dt .* (m.v[:, end] .+ cny) .* (m.v[:, end] .- m.v[:, end - 1]) ./ dy
    end

    # update! m.u, m.v
    m.u = u
    m.v = v

    # end of step t
    m.t += 1
    m
end

LinearAlgebra.adjoint(lu::ILU0Precon) = Adjoint(lu)

function LinearAlgebra.ldiv!(x::AbstractVector{M}, LU::Adjoint{M, ILU0Precon{T,N,M}}, b::AbstractVector{M}) where {T,N <: Integer,M}
    (length(b) == LU.parent.n) || throw(DimensionMismatch())
    (; n, l_colptr, l_rowval, l_nzval, u_colptr, u_rowval, u_nzval) = LU.parent
    # forward
    @inbounds for i = 1:n
        x[i] = b[i]
        for j = u_colptr[i]:u_colptr[i + 1] - 2
            x[i] -= u_nzval[j] / u_nzval[u_colptr[u_rowval[j] + 1] - 1] * x[u_rowval[j]]
        end
    end
    # backward
    @inbounds for i = n:-1:1
        x[i] /= u_nzval[u_colptr[i + 1] - 1]
        for j = l_colptr[i]:l_colptr[i + 1] - 1
            x[i] -= l_nzval[j] * x[l_rowval[j]]
        end
    end
    x
end
