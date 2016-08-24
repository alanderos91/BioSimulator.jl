type SAL <: TauLeapMethod
    # parameters
    end_time :: Float64
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

    function SAL(tf::AbstractFloat, ϵ::AbstractFloat, δ::AbstractFloat, α::AbstractFloat)
        new(tf, ϵ, δ, α,
            0.0, Float64[], Float64[], Int[], 0, 0, 0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            DEFAULT_TAULEAP)
    end
end

set_time!(algorithm::SAL, τ::AbstractFloat) = (algorithm.t = algorithm.t + τ)

##### accessors #####
epsilon(x::SAL) = x.ϵ
delta(x::SAL)   = x.δ
alpha(x::SAL)   = x.α
get_derivatives(x::SAL) = x.dxdt, x.drdt

function init!(x::SAL, Xt, r)
    c = length(Xt)
    d = size(stoichiometry(r), 2)

    setfield!(x, :dxdt,   zeros(Float64, c))
    setfield!(x, :drdt,   zeros(Float64, d))
    setfield!(x, :events, zeros(Int,     d))

    return nothing
end

function reset!(x::SAL, a::PVec)
    setfield!(x, :t, 0.0)
    setfield!(x, :ssa_steps, 0)
    setfield!(x, :leap_steps, 0)
    setfield!(x, :neg_excursions, 0)

    return nothing
end

# for use at the end of a realization
function compute_statistics!(x::SAL, i::Integer)
    n = i - 1
    setfield!(x, :avg_nsteps, cumavg(avg_nsteps(x), steps(x),     n))
    setfield!(x, :avg_nssa,    cumavg(avg_nssa(x),    ssa_steps(x), n))
    setfield!(x, :avg_nleaps,    cumavg(avg_nleaps(x),    leap_steps(x), n))
    setfield!(x, :avg_neg_excursions, cumavg(avg_neg_excursions(x), neg_excursions(x), n))
end

function step!(algorithm::SAL, Xt, r)
    a = propensities(r)

    if intensity(a) > 0
        if intensity(a) < delta(algorithm)
            τ = rand(Exponential(1 / intensity(a)))

            #setfield!(x, :ssa_steps, ssa_steps(x) + 1)
            #setfield!(x, :avg_ssa_step,  cumavg(avg_ssa_step(x),  τ, ssa_steps(x)))
            #setfield!(x, :t, time(x) + τ)
            set_time!(algorithm, τ)

            if !done(algorithm)
                μ = select_reaction(a)
                fire_reaction!(Xt, r, μ)
                update_dependent_propensities!(r, Xt, μ)
            end
        else
            τ = sal_update!(algorithm, Xt, r)
            update_all_propensities!(r, Xt)
            # setfield!(x, :leap_steps, leap_steps(x) + 1)
            # setfield!(x, :avg_leap_step,  cumavg(avg_leap_step(x),  τ, leap_steps(x)))
            # setfield!(x, :t, time(x) + τ)
            set_time!(algorithm, τ)
        end
    elseif intensity(a) == 0
      algorithm.t = algorithm.end_time
    else
      println("t = ", get_time(algorithm))
      println("a = ", a)
      println("Xt = ", Xt)
      error("intensity = ", intensity(a))
    end

    #compute_statistics!(x, τ)
    return nothing
end

function sal_update!(algorithm, Xt, r)
    dxdt, drdt = get_derivatives(algorithm)
    ϵ = epsilon(algorithm)
    α = alpha(algorithm)
    events = algorithm.events

    mean_derivatives!(dxdt, r)
    time_derivatives!(drdt, Xt, r, dxdt)

    τ = tau_leap(r, drdt, ϵ)
    τ = min(τ, end_time(algorithm) - get_time(algorithm))

    generate_events!(events, r, τ, drdt)

    while is_badleap(Xt, r, events)
        contract!(events, α)
        τ = τ * α
        #setfield!(x, :neg_excursions, neg_excursions(x) + 1)
    end

    fire_reactions!(Xt, r, events)

    return τ
end

function mean_derivatives!(dxdt, r::DenseReactionSystem)
    V = stoichiometry(r)
    a = propensities(r)

    for k in eachindex(dxdt)
        dxdt[k] = 0.0
        for j in eachindex(a)
            dxdt[k] = dxdt[k] + a[j] * V[k, j]
        end
    end

    return dxdt
end

function mean_derivatives!(dxdt, r::SparseReactionSystem)
    V = stoichiometry(r)
    a = propensities(r)

    Vj = nonzeros(V)
    ix = rowvals(V)
    c  = size(V, 2)

    fill!(dxdt, 0.0)

    for j in 1:c
        for k in nzrange(V, j)
            i = ix[k]
            dxdt[i] = dxdt[i] + a[j] * Vj[k]
        end
    end

    return dxdt
end

function time_derivatives!(drdt, Xt, r::DenseReactionSystem, dxdt)
    for i in eachindex(drdt)
        drdt[i]  = 0.0
        for k in eachindex(Xt)
            ∂r∂x_k = compute_mass_action_deriv(Xt, r, i, k)
            drdt[i] = drdt[i] + ∂r∂x_k * dxdt[k]
        end
    end
    return drdt
end

function time_derivatives!(drdt, Xt, r::SparseReactionSystem, dxdt)
    U  = coefficients(r)
    ix = rowvals(U)

    for i in eachindex(drdt)
        drdt[i]  = 0.0
        for k in nzrange(U, i)
            ∂r∂x_k = compute_mass_action_deriv(Xt, r, i, ix[k])
            drdt[i] = drdt[i] + ∂r∂x_k * dxdt[ix[k]]
        end
    end

    return drdt
end

function tau_leap(r::AbstractReactionSystem, drdt::Vector, ϵ::AbstractFloat)
    τ = Inf
    a = propensities(r)
    k = scaled_rates(r) # should we be using the unscaled rate instead?

    for j in eachindex(a)
        A = ϵ * max(a[j], k[j])
        B = abs(drdt[j])

        τ = min(τ, A / B)
    end

    return τ
end

function generate_events!(events, r, τ, drdt)
    a = propensities(r)

    for j in eachindex(a)
        λ = τ * a[j] + 0.5 * τ * τ * drdt[j]

        events[j] = rand(Poisson(max(λ, 0)))
    end

    return nothing
end

function is_badleap(Xt, r::DenseReactionSystem, events)
    V = stoichiometry(r)

    for i in eachindex(Xt)
        xi = Xt[i]
        for j in eachindex(events)
            xi = xi + events[j] * V[i, j]
            if xi < 0 return true end
        end
    end

    return false
end

function is_badleap(Xt, r::SparseReactionSystem, events)
    V  = stoichiometry(r)
    Vj = nonzeros(V)
    ix = rowvals(V)
    c  = size(V, 2)

    for j in 1:c
        xi = 0
        for k in nzrange(V, j)
            xi = Xt[ix[k]] + events[j] * Vj[k]

            if xi < 0 return true end
        end
    end

    return false
end

function contract!(events, α)
    for i in eachindex(events)
        k = 0
        for j in 1:events[i]
            k = rand() < α ? k + 1 : k
        end
        events[i] = k
    end

    return nothing
end
