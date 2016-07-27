type SAL <: TauLeapMethod
    # parameters
    T :: Float64
    ϵ :: Float64
    δ :: Float64
    α :: Float64

    # state variables
    t              :: Float64
    dxdt           :: Vector{Float64}
    drdt           :: Vector{Float64}
    events         :: Vector{Int}
    ssa_steps      :: Int
    leap_steps     :: Int
    neg_excursions :: Int

    # statistics
    avg_nsteps         :: Float64
    avg_step_size      :: Float64
    avg_nssa           :: Float64
    avg_nleaps         :: Float64
    avg_ssa_step       :: Float64
    avg_leap_step      :: Float64
    avg_neg_excursions :: Float64

    # metadata tags
    tags :: Vector{Symbol}

    function SAL(T, ϵ, δ, α)
        new(T, ϵ, δ, α,
            0.0, Float64[], Float64[], Int[], 0, 0, 0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            DEFAULT_TAULEAP)
    end
end

##### accessors #####
epsilon(x::SAL) = x.ϵ
delta(x::SAL)   = x.δ
alpha(x::SAL)   = x.α
get_derivatives(x::SAL) = x.dxdt, x.drdt

function sal(T; ϵ=0.125, δ=100.0, α=0.75, na...)
    return SAL(T, ϵ, δ, α)
end

function initialize!(x::SAL, m::Model)
    c, d = size(m)

    setfield!(x, :dxdt,   zeros(Float64, c))
    setfield!(x, :drdt,   zeros(Float64, d))
    setfield!(x, :events, zeros(Int,     d))
    return;
end

function reset!(x::SAL, m::Model)
    setfield!(x, :t, 0.0)
    setfield!(x, :ssa_steps, 0)
    setfield!(x, :leap_steps, 0)
    setfield!(x, :neg_excursions, 0)
    return;
end

# for use at the end of a realization
function compute_statistics!(x::SAL, i::Integer)
    n = i - 1
    setfield!(x, :avg_nsteps, cumavg(avg_nsteps(x), steps(x),     n))
    setfield!(x, :avg_nssa,    cumavg(avg_nssa(x),    ssa_steps(x), n))
    setfield!(x, :avg_nleaps,    cumavg(avg_nleaps(x),    leap_steps(x), n))
    setfield!(x, :avg_neg_excursions, cumavg(avg_neg_excursions(x), neg_excursions(x), n))
end

function step!(x::SAL, Xt, rs, p)
    a0 = compute_propensities!(rs, Xt, p)

    if a0 < delta(x)
        τ = rand(Exponential(1 / a0))

        setfield!(x, :ssa_steps, ssa_steps(x) + 1)
        setfield!(x, :avg_ssa_step,  cumavg(avg_ssa_step(x),  τ, ssa_steps(x)))
        setfield!(x, :t, time(x) + τ)

        if time(x) < end_time(x) && a0 > 0
            μ = select_reaction(rs, a0)
            fire_reaction!(Xt, rs, μ)
        end
    else
        τ = sal_update!(x, Xt, rs, p)
        setfield!(x, :leap_steps, leap_steps(x) + 1)
        setfield!(x, :avg_leap_step,  cumavg(avg_leap_step(x),  τ, leap_steps(x)))
        setfield!(x, :t, time(x) + τ)
    end
    compute_statistics!(x, τ)
end

function sal_update!(x, Xt, rs, p)
    dxdt, drdt = get_derivatives(x)
    ϵ = epsilon(x)
    α = alpha(x)
    events = x.events

    mean_derivatives!(dxdt, rs)
    time_derivatives!(drdt, Xt, rs, p, dxdt)

    τ = tau_leap(rs, p, drdt, ϵ)
    τ = min(τ, end_time(x) - time(x))

    generate_events!(events, rs, τ, drdt)

    while is_badleap(Xt, rs, events)
        contract!(events, α)
        τ = τ * α
        setfield!(x, :neg_excursions, neg_excursions(x) + 1)
    end

    fire_reactions!(Xt, rs, events)

    return τ
end

function mean_derivatives!(dxdt, rs::DenseReactionSystem)
    v = increments(rs)
    a = propensities(rs)
    @inbounds for k in eachindex(dxdt)
        dxdt[k] = 0.0
        @inbounds for j in eachindex(a)
            vj = v[j]
            dxdt[k] = dxdt[k] + a[j] * vj[k]
        end
    end
    return dxdt
end

function mean_derivatives!(dxdt, rs::SparseReactionSystem)
    v = increments(rs)
    a = propensities(rs)

    vj = nonzeros(v)
    idxs = rowvals(v)
    c = size(v, 2)

    fill!(dxdt, 0.0)

    @inbounds for j in 1:c
        @inbounds for k in nzrange(v, j)
            i = idxs[k]
            dxdt[i] = dxdt[i] + a[j] * vj[k]
        end
    end
    return dxdt
end

function time_derivatives!(drdt, Xt, rs::DenseReactionSystem, p, dxdt)
    @inbounds for i in eachindex(drdt)
        drdt[i]  = 0.0
        @inbounds for k in eachindex(Xt)
            ∂r∂x_k = mass_action_deriv(Xt, rs, p, i, k)
            drdt[i] = drdt[i] + ∂r∂x_k * dxdt[k]
        end
    end
    return drdt
end

function time_derivatives!(drdt, Xt, rs::SparseReactionSystem, p, dxdt)
    u  = reactants(rs)
    rv = rowvals(u)
    nz = nonzeros(u)

    @inbounds for i in eachindex(drdt)
        drdt[i]  = 0.0
        @inbounds for k in nzrange(u, i)
            ∂r∂x_k = mass_action_deriv(Xt, rs, p, i, rv[k])
            drdt[i] = drdt[i] + ∂r∂x_k * dxdt[rv[k]]
        end
    end
    return drdt
end

function tau_leap(rs, p, drdt, ϵ)
    τ = Inf
    a = propensities(rs)
    k = rates(rs)
    @inbounds for j in eachindex(a)
        kj = p[k[j]].value
        r  = a[j]

        A = ϵ * max(r, kj)
        B = abs(drdt[j])

        τ = min(τ, A / B)
    end

    return τ
end

function generate_events!(events, rs, τ, drdt)
    r = propensities(rs)
    @inbounds for j in eachindex(r)
        λ = τ * r[j] + 0.5 * τ * τ * drdt[j]
        events[j] = rand(Poisson(λ))
    end
end

function is_badleap(Xt, rs::DenseReactionSystem, events)
    v = increments(rs)

    @inbounds for i in eachindex(Xt)
        xi = Xt[i]
        @inbounds for j in eachindex(events)
            vj = v[j]
            xi = xi + events[j] * vj[i]

            if xi < 0 return true end
        end
    end

    return false
end

function is_badleap(Xt, rs::SparseReactionSystem, events)
    v = increments(rs)
    vj = nonzeros(v)
    idxs = rowvals(v)
    c = size(v, 2)

    @inbounds for j in 1:c
        xi = 0
        @inbounds for k in nzrange(v, j)
            xi = Xt[idxs[k]] + events[j] * vj[k]

            if xi < 0 return true end
        end
    end
    return false
end

function contract!(events, α)
    @inbounds for i in eachindex(events)
        k = 0
        @inbounds for j in 1:events[i]
            k = rand() < α ? k + 1 : k
        end
        events[i] = k
    end
end
