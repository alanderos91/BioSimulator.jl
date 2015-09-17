export PopulationTrace, SimulationResult, plot_trajectory

immutable PopulationState
    name::ASCIIString
    time::Float64
    value::Int
end

function Base.show(io::IO, s::PopulationState)
    @printf io " %8s   %4.4e  %9.0d\n" string(s.name) s.time s.value
end



immutable PopulationTrace
    states::Vector{PopulationState}
end

function PopulationTrace()
    return PopulationTrace(PopulationState[])
end

function PopulationTrace(x::Species)
    state = PopulationState(x.id, 0.0, x.pop)
    return PopulationTrace([state])
end

Base.push!(t::PopulationTrace, s::PopulationState) = push!(t.states, s)
Base.getindex(t::PopulationTrace, i::Integer) = getindex(t.states, i)

function Base.setindex!(t::PopulationTrace, s::PopulationState, i::Integer)
    setindex!(t.states, s, i)
end

function Base.isempty(t::PopulationTrace)
    return isempty(t.states)
end

function Base.show(io::IO, t::PopulationTrace)
    @printf io "   Name        Time       Value  \n"
    @printf io "----------  ----------  ---------\n"
    for state in t.states
        show(io, state)
    end
end

function update_traces!(traces::Dict{ASCIIString, PopulationTrace}, t::Float64, spcs::Vector{Species}, store_trace::Bool)
    for s in spcs
        s.istracked ? update!(traces[s.id], t, s, store_trace) : Nothing
    end
end

function update!(tr, t, s, store_trace)
    state = PopulationState(s.id, t, s.pop)

    if store_trace
        push!(tr, state)
    end
end

function plot_trajectory(tr::PopulationTrace)
    if !isempty(tr)
        name = tr.states[1].name
        tval = Float64[]
        pval = Int[]

        for state in tr.states
            push!(tval, state.time)
            push!(pval, state.value)
        end

        return Gadfly.plot(x=tval, y=pval, Geom.line, Guide.xlabel("Time"), Guide.ylabel("Population"), Guide.title(string(name, " vs Time")))
    else

    end
end



immutable SimulationResult
    algorithm::ASCIIString
    traces::Dict{ASCIIString, PopulationTrace}
    metadata::Dict{}
end
