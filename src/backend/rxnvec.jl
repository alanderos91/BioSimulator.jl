##### Reaction Channel #####

immutable ReactionChannel
    id::Symbol
    rate::Symbol
    pre::Vector{Int}
    post::Vector{Int}

    function ReactionChannel(id, rate, pre, post)
        if any(pre .< 0) || any(post .< 0)
          error("stoichiometric coefficients must be positive.")
        end
        return new(id, rate, pre, post)
    end
end

function ReactionChannel(r::Reaction, id2ind)
    pre  = Array(Int, length(id2ind))
    post = Array(Int, length(id2ind))

    for (id, i) in id2ind
        pre[i]  = get(r.reactants, id, 0)
        post[i] = get(r.products,  id, 0)
    end

    return ReactionChannel(r.id, r.rate, pre, post)
end

increment(rc::ReactionChannel, i) = rc.post[i] - rc.pre[i]
rate(rc::ReactionChannel) = rc.rate

##### Reaction Vector #####

type ReactionVector
    reactions::Vector{ReactionChannel}
    intensity::Float64
    propensities::Vector{Float64}
end

function ReactionVector(rxns::Vector{ReactionChannel})
  return ReactionVector(rxns, Inf, zeros(Float64, length(rxns)))
end

##### accessors #####
reactions(rv::ReactionVector)    = rv.reactions
intensity(rv::ReactionVector)    = rv.intensity
propensities(rv::ReactionVector) = rv.propensities
reaction(rv::ReactionVector, i)  = rv.reactions[i]

##### iteration #####
start(rv::ReactionVector)       = start(rv.propensities)
done(rv::ReactionVector, state) = done(rv.propensities, state)
next(rv::ReactionVector, state) = next(rv.propensities, state)
enumerate(rv::ReactionVector)   = enumerate(rv.propensities)
eachindex(rv::ReactionVector)   = eachindex(rv.propensities)

##### general collection #####
isempty(rv::ReactionVector) = isempty(rv.propensities)
empty!(rv::ReactionVector)  = empty!(rv.propensities)
length(rv::ReactionVector)  = length(rv.propensities)
endof(rv::ReactionVector)   = endof(rv.propensities)

##### indexing #####
getindex(rv::ReactionVector, key...)         = getindex(rv.propensities, key...)
setindex!(rv::ReactionVector, value, key...) = setindex!(rv.propensities, value, key...)

##### sorting #####
function sort!(rv::ReactionVector;
               alg::Base.Algorithm=defalg(v),
               lt=isless,
               by=identity,
               rev::Bool=false,
               order::Base.Ordering=Base.Forward)
    rxns = reactions(rv)
    prop = propensities(rv)
    ix   = sortperm(prop, alg=alg, lt=lt, by=by, order=order)

    for i in eachindex(ix)
        k = ix[i]
        r = rxns[i]
        p = prop[i]

        rxns[i] = rxns[k]
        rxns[k] = r

        prop[i] = prop[k]
        prop[k] = p
    end
    return rv
end

##### mass-action kinetics ######
function mass_action(r::ReactionChannel, x::Vector{Int}, params::Parameters)
  c = params[r.rate]
  return mass_action(c, r.pre, x)
end

function mass_action_deriv(r::ReactionChannel, x::Vector{Int}, params::Parameters, k::Int)
  c = params[r.rate]
  return mass_action_deriv(c, r.pre, x, k)
end

function intensity(r::ReactionChannel, x::Vector{Int}, param::Parameters)
    mass_action(r, x, param)
end

function compute_intensities!(rv::ReactionVector, x::Vector{Int}, param::Parameters)
  a0  = 0.0
  rxn = reactions(rv)
  for j in eachindex(rxn)
    aj = intensity(rxn[j], x, param)
    rv[j] = aj
    a0 = a0 + aj
  end
  rv.intensity = a0
end

##### misc #####
fill!(rv::ReactionVector, value) = fill!(rv.propensities, value)
