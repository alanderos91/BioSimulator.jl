using IteratorInterfaceExtension, TableTraits

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
function IteratorInterfaceExtensions.getiterator(xw::SamplePath)
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

# use generated function to help with type inference?
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
