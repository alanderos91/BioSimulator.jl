using IteratorInterfaceExtensions, TableTraits

import IteratorInterfaceExtensions: isiterable, getiterator
import TableTraits: isiterabletable
import Base: length, eltype, iterate

# T stores the type info for a NamedTuple
# S stores the SamplePath type
struct SamplePathIterator{T,S}
    xw::S
end

# required methods for SamplePath
IteratorInterfaceExtensions.isiterable(xw::SamplePath) = true
TableTraits.isiterabletable(xw::SamplePath) = true

# TODO: how should we handle Configuration data?
function IteratorInterfaceExtensions.getiterator(xw::SamplePath{T}) where T <: Number
    t_type = eltype(xw.t)

    # format:
    # t X1 X2 X3 ...
    value_types = Type[]
    value_names = Symbol[]

    push!(value_types, t_type)
    push!(value_names, :t)

    # use the first record as a representative
    example = xw.u[1]

    # fill in the types and names for each component; should be the same
    # should be the same (for now)
    for k in 1:length(example)
        push!(value_types, eltype(example))
        push!(value_names, Symbol("X$(k)"))
    end

    # construct our iterator
    iter = SamplePathIterator{NamedTuple{(value_names...,), Tuple{value_types...}}, typeof(xw)}(xw)

    return iter
end

# iteration interface implementation for SamplePathIterator
Base.length(iter::SamplePathIterator) = length(iter.xw)
Base.eltype(iter::SamplePathIterator{T,S}) where {T,S} = T

@generated function Base.iterate(iter::SamplePathIterator{T,S}, state=1) where {T,S}
    record = []

    # add t value for the record
    push!(record, :(iter.xw.t[state]))

    # populate the record with values for each component
    n = length(fieldnames(T))
    for i in 1:n-1
        push!(record, :(iter.xw.u[state][$i]))
    end

    # function body
    # each iterate is a named tuple
    quote
        if state > length(iter)
            return nothing
        else
            return $(T)(($(record...),)), state+1
        end
    end
end

# workaround for broken integrations in IterableTables.jl
tablefy(xw::SamplePath) = IteratorInterfaceExtensions.getiterator(xw)

# T stores the type info for a NamedTuple
# S stores the SamplePath type
struct EnsembleIterator{T,S}
    x::S
    trials::Int
end

# required methods for SamplePath
IteratorInterfaceExtensions.isiterable(x::Ensemble) = true
TableTraits.isiterabletable(x::Ensemble) = true

# TODO: how should we handle Configuration data?
function IteratorInterfaceExtensions.getiterator(x::Ensemble{T}) where T <: Number
    t_type = eltype(x[1].t)

    # format:
    # trial t X1 X2 X3 ...
    value_types = Type[]
    value_names = Symbol[]

    # first column is iteration number
    push!(value_types, Int)
    push!(value_names, :trial)

    # second column is a timestamp
    push!(value_types, t_type)
    push!(value_names, :t)

    # use the first record as a representative
    example = x[1].u[1]

    # fill in the types and names for each component; should be the same
    # should be the same (for now)
    for k in 1:length(example)
        push!(value_types, eltype(example))
        push!(value_names, Symbol("X$(k)"))
    end

    # construct our iterator
    iter = EnsembleIterator{NamedTuple{(value_names...,), Tuple{value_types...}}, typeof(x)}(x, length(x))

    return iter
end

# iteration interface implementation for SamplePathIterator
Base.length(iter::EnsembleIterator) = sum(length, iter.x)
Base.eltype(iter::EnsembleIterator{T,S}) where {T,S} = T

@generated function Base.iterate(iter::EnsembleIterator{T,S}, state=(1,1)) where {T,S}
    record = []
    # add trial number
    push!(record, :(state[1]))

    # add t value for the record
    push!(record, :(iter.x[state[1]].t[state[2]]))

    # populate the record with values for each component
    n = length(fieldnames(T))

    for i in 1:n-2
        push!(record, :(iter.x[state[1]].u[state[2]][$i]))
    end

    # function body
    # each iterate is a named tuple
    quote
        if state[1] > iter.trials
            return nothing
        else
            # get number of states for current trajectory
            m = length(iter.x[state[1]])

            # if we've touched all the values for this trajectory,
            # move on to the next one
            if state[2] ≥ m
                next_state = (state[1]+1, 1)
            else
            # otherwise, advance to the next time point
                next_state = (state[1], state[2]+1)
            end

            return $(T)(($(record...),)), next_state
        end
    end
end

# workaround for broken integrations in IterableTables.jl
tablefy(x::Ensemble) = IteratorInterfaceExtensions.getiterator(x)
