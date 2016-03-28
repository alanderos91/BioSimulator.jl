type SAL <: Algorithm
    # parameters
    T::Float64
    ϵ::Float64
    δ::Float64
    α::Float64

    # state variables
    t::Float64
    dxdt::Vector{Float64}
    drdt::Vector{Float64}
    events::Vector{Int}
    ssa_steps::Int
    sal_steps::Int

    # statistics
    avg_nsteps::Float64
    avg_step_size::Float64
    avg_ssa::Float64
    avg_sal::Float64
    avg_ssa_step::Float64
    avg_sal_step::Float64

    # metadata tags
    tags::Vector{Symbol}

    function SAL(T, ϵ, δ, α)
        new(T, ϵ, δ, α,
            0.0, Float64[], Float64[], Int[], 0, 0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            [:avg_nsteps, :avg_step_size, :avg_nsteps, :avg_step_size, :avg_ssa, :avg_sal, :avg_ssa_step, :avg_sal_step])
    end
end

##### accessors #####
epsilon(x::SAL) = x.ϵ
delta(x::SAL)   = x.δ
alpha(x::SAL)   = x.α

get_derivatives(x::SAL) = x.dxdt, x.drdt
events(x::SAL)          = x.events
steps(x::SAL)           = x.ssa_steps + x.sal_steps
ssa_steps(x::SAL)       = x.ssa_steps
sal_steps(x::SAL)       = x.sal_steps

avg_nsteps(x::SAL)   = x.avg_nsteps
avg_ssa(x::SAL)      = x.avg_ssa
avg_sal(x::SAL)      = x.avg_sal
avg_ssa_step(x::SAL) = x.avg_ssa_step
avg_sal_step(x::SAL) = x.avg_sal_step

function sal(T; ϵ=0.125, δ=100.0, α=0.75, na...)
    return SAL(T, ϵ, δ, α)
end

function initialize!(x::SAL, m::Model)
    c = length(species(m))
    d = length(reactions(m))

    setfield!(x, :dxdt,   zeros(Float64, c))
    setfield!(x, :drdt,   zeros(Float64, d))
    setfield!(x, :events, zeros(Int,     d))
end

function reset!(x::SAL, m::Model)
    setfield!(x, :t, 0.0)
    setfield!(x, :ssa_steps, 0)
    setfield!(x, :sal_steps, 0)
    return;
end

function compute_statistics!(x::SAL, τ::Float64)
    setfield!(x, :avg_step_size, cumavg(avg_step_size(x), τ, steps(x)))
end

# for use at the end of a realization
function compute_statistics!(x::SAL, i::Integer)
    setfield!(x, :avg_nsteps, cumavg(avg_nsteps(x), steps(x),     i))
    setfield!(x, :avg_ssa,    cumavg(avg_ssa(x),    ssa_steps(x), i))
    setfield!(x, :avg_sal,    cumavg(avg_sal(x),    sal_steps(x), i))
end

function call(x::SAL, m::Model)
    a0 = compute_propensities!(m)

    if a0 < delta(x)
        τ = ssa_update!(m, x)
        setfield!(x, :ssa_steps, ssa_steps(x) + 1)
        setfield!(x, :avg_ssa_step,  cumavg(avg_ssa_step(x),  τ, ssa_steps(x)))
    else
        τ = sal_update!(m, x)
        setfield!(x, :sal_steps, sal_steps(x) + 1)
        setfield!(x, :avg_sal_step,  cumavg(avg_sal_step(x),  τ, sal_steps(x)))
    end

    setfield!(x, :t, time(x) + τ)
    compute_statistics!(x, τ)
end

function sal_update!(m::Model, x)
    dxdt, drdt = get_derivatives(x)
    ϵ = epsilon(x)
    α = alpha(x)
    events = x.events

    mean_derivatives!(dxdt, m)
    time_derivatives!(drdt, m, dxdt)

    τ = tau_leap(m, drdt, ϵ)
    τ = min(τ, end_time(x) - time(x))

    generate_events!(events, m, τ, drdt)

    while is_badleap(m, events)
        contract!(events, α)
        τ = τ * α
    end

    fire_reactions!(m, events)

    return τ
end

function mean_derivatives!(dxdt, m)
    rxns = reactions(m)

    for k in eachindex(dxdt)
        dxdt[k] = 0.0
        for j in eachindex(rxns)
            r = reaction(rxns, j)
            dxdt[k] = dxdt[k] + rxns[j] * increment(r, k)
        end
    end
    return dxdt
end

function time_derivatives!(drdt, m, dxdt)
    Xt   = species(m)
    rxns = reactions(m)
    p    = parameters(m)

    for i in eachindex(drdt)
        drdt[i]  = 0.0
        rxn = reaction(rxns, i)
        for k in eachindex(Xt)
            ∂r∂x_k = mass_action_deriv(rxn, Xt, p, k)
            drdt[i] = drdt[i] + ∂r∂x_k * dxdt[k]
        end
    end
    return drdt
end

function tau_leap(m, drdt, ϵ)
    τ    = Inf
    rxns = reactions(m)
    p    = parameters(m)
    for j in eachindex(rxns)
        r   = rxns[j]
        key = rate(reaction(rxns, j))
        k   = p[key]

        a = ϵ * max(r, k)
        b = abs(drdt[j])

        τ = min(τ, a / b)
    end

    return τ
end

function generate_events!(events, m, τ, drdt)
    rxn = reactions(m)
    for j in eachindex(rxn)
        λ = τ * rxn[j] + 0.5 * τ * τ * drdt[j]
        events[j] = rand(Poisson(λ))
    end
end

function is_badleap(m, events)
    Xt   = species(m)
    rxns = reactions(m)

    for i in eachindex(Xt)
        xi = Xt[i]
        for j in eachindex(rxns)
            xi = xi + events[j] * increment(reaction(rxns, j), i)

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
end
