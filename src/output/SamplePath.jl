struct SamplePath{T,N,A,B} <: AbstractDiffEqArray{T, N}
    u::A # A <: AbstractVector{<: AbstractArray{T, N - 1}}
    t::B
end

Ensemble{T,N,A,B} = Vector{SamplePath{T,N,A,B}}

SamplePath(xs::AbstractVector{T}, ts, dims::NTuple{N}) where {T, N} = SamplePath{eltype(T), N, typeof(xs), typeof(ts)}(xs, ts)
# Assume that the first element is representative all all other elements
SamplePath(xs::AbstractVector, ts::AbstractVector) = SamplePath(xs, ts, (size(xs[1])..., length(xs)))

# @inline Base.getindex(xw::SamplePath{T, N}, I::AbstractArray{Int}) where {T, N} = SamplePath(xw.u[I], xw.t[I])

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

# SamplePath() = SamplePath{Float64,Int}()

# Ensemble(ntrials) = [SamplePath() for i in 1:ntrials]

# function Base.sizehint!(xw::SamplePath, n)
#     sizehint!(xw.t, n)
#     sizehint!(xw.u, n)
#
#     return xw
# end

function update!(xw::SamplePath, t, x, save_points)
    j = searchsortedlast(save_points, t)
    i = j

    # backtrack in case the last event jumped over save points
    while (i > 1) && !(save_points[i] in xw.t)
        i -= 1
    end

    # fill in the data for each save point
    for k in i+1:j
        push!(xw.u, copy(x))
        push!(xw.t, save_points[k])
    end

    return xw
end

function update!(xw::SamplePath, t, x, save_points::Nothing)
    push!(xw.u, copy(x))
    push!(xw.t, t)

    return xw
end

##### interpolation
function get_regular_path(xw::SamplePath, tfinal, epochs) where {T1,T2}
    max_epoch = epochs + 1

    if xw.u[1] isa Number
        x0 = [xw.u[1]]
    else
        x0 = xw.u[1]
    end

    ts = collect(range(0.0, stop = tfinal, length = max_epoch))
    new_xw = SamplePath([x0], ts)
    epoch = 2

    # copy data from the sample path
    for j in 2:length(xw)
        @inbounds t = xw.t[j]
        @inbounds x = xw.u[j]

        while (epoch <= max_epoch) && (t >= ts[epoch])
            push!(new_xw.u, x)
            epoch = epoch + 1
        end
    end

    # fill in the regular path
    while epoch <= max_epoch
        x = xw.u[end]
        push!(new_xw.u, x)
        epoch = epoch + 1
    end

    return new_xw
end

function get_regular_ensemble(x::Ensemble, tfinal, epochs)
    return [get_regular_path(x[i], tfinal, epochs) for i in eachindex(x)]
end

##### plotting recipes #####

@recipe function f(xw::SamplePath)
    seriestype --> :steppre

    xw.t, xw'
end

@recipe function f(xw::SamplePath{T,1}) where {T}
    seriestype --> :steppre

    xw.t, xw.u
end

@recipe function f(ens::Ensemble;
    summary = :mean,
    epochs = 100,
    error_style = :bars,
    timepoint = 1,
    vars = nothing)
    # regularize the sample paths
    tfinal = ens[1].t[end]

    reg = get_regular_ensemble(ens, tfinal, epochs)

    if summary == :mean
        # extract the series data
        ts = reg[1].t
        xs = transpose(mean(reg))
        xstd = transpose(std(reg))

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

        seriestype --> :scatter
        marker_z --> reg[1].t

        xs, ys
    else
        error("Unknown option for summary_type: $(summary_type)")
    end
end
