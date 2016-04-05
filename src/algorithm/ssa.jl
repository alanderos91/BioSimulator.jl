"""
```
SSA(T)
```

Gillespie's Stochastic Simulation Algorithm. Simulates a system of coupled species and reactions by simulating each reaction event. The algorithm computes the time to the next reaction as a random exponential deviate with mean equal to the cumulative reaction intensity. It then scales the cumulative intensity by a uniform random number and conducts a linear search on the reaction propensities to select a reaction.

### Arguments
- `T`: The simulation end time.
"""
type SSA <: ExactMethod
    # parameters
    T :: Float64

    # state variables
    t     :: Float64
    steps :: Int

    # statistics
    avg_nsteps    :: Float64
    avg_step_size :: Float64

    # metadata tags
    tags :: Vector{Symbol}

    function SSA(T)
        new(T, 0.0, 0, 0.0, 0.0, DEFAULT_EXACT)
    end
end

##### constructor wrapper #####
function ssa(T; na...)
    return SSA(T)
end

function step!(x::SSA, Xt, rs, p)
    a0 = compute_propensities!(rs, Xt, p)
    τ  = rand(Exponential(1 / a0))

    # update algorithm variables
    setfield!(x, :t,     time(x) + τ)
    setfield!(x, :steps, steps(x) + 1)
    compute_statistics!(x, τ)

    if time(x) < end_time(x) && a0 > 0
        μ = select_reaction(rs, a0)
        fire_reaction!(Xt, rs, μ)
    end
    return;
end


function select_reaction(rs::AbstractReactionSystem, a0)
    jump = a0 * rand()
    a = propensities(rs)
    b = 0.0
    @inbounds for j in eachindex(a)
        b = b + a[j]
        if b >= jump return j end
    end
    return 0
end
