"""
```
FRM(T)
```

Gillespie's First Reaction Method. Statistically equivalent to `SSA`, but slower computationally. The algorithm computes the time to the next reaction as the minimum the relative firing times for each reaction.

### Arguments
- `T`: The simulation end time.
"""
type FRM <: ExactMethod
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

    function FRM(T)
        new(T, 0.0, 0, 0.0, 0.0, DEFAULT_EXACT)
    end
end

##### constructor wrapper #####
function frm(T; na...)
    return FRM(T)
end

function step!(x::FRM, Xt, rs, p)
    a0 = compute_propensities!(rs, Xt, p)
    τ, μ = select_reaction(rs)

    # update algorithm variables
    setfield!(x, :t,     time(x) + τ)
    setfield!(x, :steps, steps(x) + 1)

    if time(x) < end_time(x) && a0 > 0
        if μ > 0
            fire_reaction!(Xt, rs, μ)
        else
            error("no reaction ocurred!")
        end
    end

    compute_statistics!(x, τ)
end

function select_reaction(rs::AbstractReactionSystem)
    τ = Inf
    μ = 0
    a = propensities(rs)
    @inbounds for j in eachindex(a)
        τj = rand(Exponential(1 / a[j]))
        if τj < τ
            τ = τj
            μ = j
        end
    end
    return τ, μ
end
