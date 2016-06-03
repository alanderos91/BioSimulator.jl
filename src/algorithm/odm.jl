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
    T          :: Float64
    init_steps :: Int
    init_iters :: Int

    # state variables
    t     :: Float64
    steps :: Int

    # statistics
    avg_nsteps    :: Float64
    avg_step_size :: Float64

    # metadata tags
    tags :: Vector{Symbol}

    function ODM(T, init_steps, init_iters)
        new(T, init_steps, init_iters, 0.0, 0, 0.0, 0.0, DEFAULT_EXACT)
    end
end

function odm(T; init_steps=100, init_iters=1, na...)
    return ODM(T, init_steps, init_iters)
end

function initialize!(x::ODM, m::Model)

    # Presimulate to sort reactions according to multiscale property. This will modify Xt and rs...
    presimulate!(x, m)
    return;
end

function reset!(x::ODM, m::Model)
    Xt = species(m)
    rs = reactions(m)
    p  = parameters(m)

    setfield!(x, :t, 0.0)
    setfield!(x, :steps, 0)
    compute_propensities!(rs, Xt, p)

    return;
end

function step!(x::ODM, Xt, rs, p)
    a0 = intensity(propensities(rs))
    τ  = rand(Exponential(1 / a0))

    # update algorithm variables
    setfield!(x, :t,     time(x) + τ)
    setfield!(x, :steps, steps(x) + 1)

    if time(x) < end_time(x) && a0 > 0
        μ = select_reaction(rs, a0)
        fire_reaction!(Xt, rs, μ)
        compute_propensities!(rs, Xt, p)
    end

    compute_statistics!(x, τ)
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
    a = propensities(rs).a
    for i in eachindex(ix)
        swapcols!(u, i, ix[i])
        swapcols!(v, i, ix[i])
        swapcols!(a, i, ix[i])
    end
    return rs
end

function swapcols!(A, i, j)
     for k = 1:size(A, 1)
           A[k,i], A[k,j] = A[k,j], A[k,i]
     end
     return A
end

function swapcols!(v::Vector, i, j)
    v[i], v[j] = v[j], v[i]
    return v
end
