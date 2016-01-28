abstract Observer
id(o::Observer)              = o.id
states(o::Observer)          = o.states
is_preallocated(o::Observer) = o.is_preallocated

immutable SpeciesObserver <: Observer
    id::Symbol
    species::Vector{Int}
    states::Vector{Int}
    is_preallocated::Bool
    index::Int

    function SpeciesObserver(id, species, index, n)
        states = zeros(Int, n)
        flag   = !isempty(states)
        new(id, species, states, flag, index)
    end
end

species(o::SpeciesObserver) = o.species
index(o::SpeciesObserver)   = o.index

function update!(observer::SpeciesObserver, j)
    i       = BioSimulator.index(observer)
    states  = BioSimulator.states(observer)
    species = BioSimulator.species(observer)

    if is_preallocated(observer)
        states[j] = species[i]
    else
        push!(states, species[i])
    end
end

immutable PropensityObserver <: Observer
    id::Symbol
    reaction::ReactionChannel
    states::Vector{Float64}
    is_preallocated::Bool

    function PropensityObserver(id, reaction, n)
        states = zeros(Float64, n)
        flag   = !isempty(states)
        new(id, reaction, states, flag)
    end
end

reaction(o::PropensityObserver) = o.reaction

function update!(observer::PropensityObserver, j)
    states   = BioSimulator.states(observer)
    reaction = BioSimulator.reaction(observer)

    if is_preallocated(observer)
        states[j] = reaction.propensity
    else
        push!(states, reaction.propensity)
    end
end

immutable TimeObserver <: Observer
    id::Symbol
    states::Vector{Float64}
    is_preallocated::Bool

    function TimeObserver(id, n)
        states = zeros(Float64, n)
        flag   = !isempty(states)
        new(id, states, flag)
    end
end

function update!(observer::TimeObserver, j, t)
    states = BioSimulator.states(observer)

    if is_preallocated(observer)
        states[j] = t
    else
        push!(states, t)
    end
end

immutable Overseer
    s_observers::Vector{SpeciesObserver}
    r_observers::Vector{PropensityObserver}
    t_observer::TimeObserver

    function Overseer(observer)
        new(SpeciesObserver[], PropensityObserver[], observer)
    end
end

function notify!(overseer, j, t)
    notify!(overseer, j)
    update!(overseer.t_observer, j, t)

    return overseer
end

function notify!(overseer, j)
    for o in overseer.s_observers
        update!(o, j)
    end

    for o in overseer.r_observers
        update!(o, j)
    end

    return overseer
end
