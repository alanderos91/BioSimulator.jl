"""
```
ODM(T)
```

Optimized Direct Method. Same as SSA, except the system is pre-simulated in order to sort reaction propensities from increasing to decreasing. This improves the search on the propensity vector when selecting the next reaction to fire.

### Arguments
- `T`: The simulation end time.

### Optional Arguments
- `init_steps`: Number of time steps to pre-simulate.
- `init_iters`: Number of iterations to simulate.
"""
type ODM <: ExactMethod
    # parameters
    end_time   :: Float64
    init_steps :: Int
    init_iters :: Int

    # state variables
    t     :: Float64
    steps :: Int

    # statistics
    avg_nsteps    :: Float64
    avg_stepsz :: Float64

    # metadata tags
    tags :: Vector{Symbol}

    function ODM(T, init_steps, init_iters)
        new(T, init_steps, init_iters, 0.0, 0, 0.0, 0.0, DEFAULT_EXACT)
    end
end

function initialize!(x::ODM, m::Model)

    presimulate!(x, m)
    return;
end

# function reset!(x::ODM, m::Model)
#     Xt = species(m)
#     rs = reactions(m)
#     p  = parameters(m)
#
#     setfield!(x, :t, 0.0)
#     setfield!(x, :steps, 0)
#     compute_propensities!(rs, Xt, p)
#
#     return;
# end

function step!(algorithm::ODM, Xt, r)
    a = propensities(r)

    τ  = select_reaction(a)

    # update algorithm variables
    set_time!(algorithm, τ)

    if !done(algorithm) && intensity(a) > 0
        μ = select_reaction(a)
        fire_reaction!(Xt, r, μ)
        update_propensities!(r, Xt, μ)
    end

    return nothing
end

function presimulate!(x, m)
    itr = x.init_iters
    n   = x.init_steps

    Xt = species(m)
    rs = reactions(m)
    p  = parameters(m)

    a = propensities(rs)
    events = zeros(Float64, length(a))

    for i = 1:itr
        reset!(m)
        for k = 1:n
            a0 = compute_propensities!(rs, Xt, p)
            μ  = select_reaction(rs, a0)
            fire_reaction!(Xt, rs, μ)
            events[μ] = events[μ] + 1
        end
    end

    for i in eachindex(a)
        @inbounds a[i] = events[i] / itr
    end

    ix = sortperm(a.a)
    _sort!(rs, ix)
end

function _sort!(rs, ix)
    u = reactants(rs)
    v = increments(rs)

    swapcols!(u, ix)
    swapcols!(v, ix)

    return rs
end

function swapcols!(A, ix)
    for i in eachindex(ix)
        j = ix[i]
         for k = 1:size(A, 1)
               A[k,i], A[k,j] = A[k,j], A[k,i]
         end
     end
     return A
end

function swapcols!(v::Vector, ix)
    v = v[ix]
    return v
end
