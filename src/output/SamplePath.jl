struct SamplePath{T,N,A,B} <: AbstractDiffEqArray{T, N}
    u::A # A <: AbstractVector{<: AbstractArray{T, N - 1}}
    t::B
end

Ensemble{T,N,A,B} = Vector{SamplePath{T,N,A,B}}

SamplePath(xs::AbstractVector{T}, ts, dims::NTuple{N}) where {T, N} = SamplePath{eltype(T), N, typeof(xs), typeof(ts)}(xs, ts)
# Assume that the first element is representative all all other elements
SamplePath(xs::AbstractVector, ts::AbstractVector) = SamplePath(xs, ts, (size(xs[1])..., length(xs)))

##### pretty printing
function Base.show(io::IO, xw::SamplePath)
    print(io,"t: "); show(io, xw.t)
    println(io);
    print(io,"x: "); show(io, xw.u)
    nothing
end

function Base.show(io::IO, m::MIME"text/plain", xw::SamplePath)
    print(io,"t: "); show(io,m,xw.t)
    println(io)
    print(io,"x: "); show(io,m,xw.u)
    nothing
end

##### update! logic
function update!(xw::SamplePath, t, state, save_points)
    j = searchsortedlast(save_points, t)
    i = j

    # backtrack in case the last event jumped over save points
    while (i > 1) && !(save_points[i] in xw.t)
        i -= 1
    end

    # fill in the data for each save point
    for k in i+1:j
        push!(xw.u, __extract(state))
        push!(xw.t, save_points[k])
    end

    return xw
end

function update!(xw::SamplePath, t, x, save_points::Nothing)
    push!(xw.u, __extract(x))
    push!(xw.t, t)

    return xw
end

##### extracting state
__extract(state::Vector{T}) where T <: Int = copy(state)
__extract(state::Lattice) = Configuration(state)

##### initializing output
function build_output(state, model)
    xw = SamplePath([__extract(state)], [0.0])
    sizehint!(xw, 1_000)

    return xw
end

##### interpolation
function get_regular_path(xw::SamplePath, save_points)
    yw = SamplePath([__extract(xw.u[1])], [save_points[1]])

    # copy data from the sample path
    for j in 2:length(xw.u)
        update!(yw, xw.t[j], xw.u[j], save_points)
    end

    return yw
end

function get_regular_ensemble(x::Ensemble, save_points)
    return [get_regular_path(x[i], save_points) for i in eachindex(x)]
end

##### plotting recipes #####
__make_indices(vars::Integer) = vars
__make_indices(vars::NTuple{N,Integer}) where N = (N > 1 ? collect(vars) : vars[1])
__make_indices(vars::AbstractVector{T}) where T <: Integer = vars

@recipe function f(xw::SamplePath{T};
    summary = :trajectory,
    vars = nothing) where T <: Number

    if summary == :trajectory
        seriestype --> :steppre
        n = length(xw[1])
        if vars == nothing
            xw.t, xw'
        else
            idxs = __make_indices(vars)

            if idxs isa Integer
                xw.t, xw[idxs,:]
            else
                xw.t, xw[idxs,:]'
            end
        end
    elseif summary == :phase
        idxs = (vars == nothing) ? (1,2) : (vars[1], vars[2])
        xs = [xw[idxs[1],t] for t in eachindex(xw.t)]
        ys = [xw[idxs[2],t] for t in eachindex(xw.t)]

        legend --> nothing

        @series begin
            seriestype --> :scatter
            marker_z --> xw.t

            xs, ys
        end

        @series begin
            seriestype --> :path
            xs, ys
        end
    end
end

@recipe function f(ens::Ensemble;
    summary = :mean,
    save_points = nothing,
    error_style = :bars,
    timepoint = 1,
    vars = nothing)

    if save_points == nothing
        m = length(ens[1].t)
        if any(xw -> length(xw) != m, ens)
            t_values = 0:ens[1].t[end]
        else
            t_values = ens[1].t
        end
    else

    end

    reg = get_regular_ensemble(ens, t_values)

    if summary == :mean
        # extract the series data
        ts = reg[1].t

        tmp_xs = transpose(mean(reg))
        tmp_xstd = transpose(std(reg))

        if vars == nothing
            xs = tmp_xs
            xstd = tmp_xstd
        else
            idxs = __make_indices(vars)

            xs = tmp_xs[:, idxs]
            xstd = tmp_xstd[:, idxs]
        end

        # make the default a scatter plot
        # and add bars to represent standard deviation
        seriestype --> :scatter

        if error_style == :bars
            yerrorbar  --> xstd
        elseif error_style == :ribbon
            ribbon --> xstd
        else
            error("Unknown option for error_style: $(error_style)")
        end

        ts, xs
    elseif summary == :histogram
        t = max(1, searchsortedlast(reg[1].t, timepoint))
        default_idxs = 1:length(reg[1].u[1])
        idxs = (vars == nothing) ? default_idxs : vars

        seriestype --> :histogram

        xs = [reg[n][k,t] for n in eachindex(reg), k in idxs]
    elseif summary == :phase
        tmp = mean(reg)
        idxs = (vars == nothing) ? (1,2) : (vars[1], vars[2])
        xs = [tmp[idxs[1],t] for t in eachindex(reg[1].t)]
        ys = [tmp[idxs[2],t] for t in eachindex(reg[1].t)]

        legend --> nothing

        @series begin
            seriestype --> :scatter
            marker_z --> reg[1].t

            xs, ys
        end

        @series begin
            seriestype --> :path
            xs, ys
        end
    else
        error("Unknown option for summary_type: $(summary_type)")
    end
end
