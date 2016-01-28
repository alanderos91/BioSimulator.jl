import Base.show

immutable SimulationOutput
    species_data::DataFrame
    propensity_data::DataFrame
end

get_species_data(so::SimulationOutput)    = so.species_data
get_propensity_data(so::SimulationOutput) = so.propensity_data

function compile_data(overseer)
    species_data    = DataFrame()
    propensity_data = DataFrame()

    t_observer  = overseer.t_observer
    s_observers = overseer.s_observers
    r_observers = overseer.r_observers

    if !isempty(t_observer.states)
        species_data[t_observer.id]    = t_observer.states
        propensity_data[t_observer.id] = t_observer.states
    end

    for o in s_observers
        species_data[o.id] = o.states
    end

    for o in r_observers
        propensity_data[o.id] = o.states
    end

    return SimulationOutput(species_data, propensity_data)
end

function Base.show(io::IO, x::SimulationOutput)
    @printf io "[species data] %d x %d\n" size(x.species_data, 1) size(x.species_data, 2)
    @printf io "[propensity data] %d x %d\n" size(x.propensity_data, 1) size(x.propensity_data, 2)
end
