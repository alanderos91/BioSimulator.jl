function extract_index_ids(id2ind, select)
    if isempty(select)
        ix  = collect(values(id2ind))
        ids = collect(keys(id2ind))
    else
        ix  = Int[ id2ind[symbol(id)] for id in select ]
        ids = map(symbol, select)
    end

    return ix, ids
end

@userplot MeanTrajectory

@recipe function plot(mt::MeanTrajectory; select=[])
    t, μ, σ, ids = meanplot(mt.args[1], select)

    # global attributes
    legend -->  true
    grid   --> false
    xguide --> "time"
    yguide --> "population mean"
    xlims  --> (t[1], t[end])
    ylims  --> (0.0, Inf)
    fillalpha --> 0.3
    seriestype --> :path
    label --> ids'
    ribbon --> σ

    t, μ
end


function meanplot(result::PartialHistory, select)
    t      = result.t
    data   = result.data
    id2ind = result.id2ind

    ix, ids = extract_index_ids(id2ind, select)

    tmp = data[ix, :, :]

    μ = transpose(reshape(mean(tmp, 3), size(tmp, 1, 2)))
    σ = transpose(reshape(std(tmp, 3), size(tmp, 1, 2)))

    return t, μ, σ, ids
end

@userplot FreqHistogram

@recipe function plot(fh::FreqHistogram; select=[])
    counts, ids = freqplot(fh.args[1], fh.args[2], select)

    # global attributes
    legend -->  true
    grid   --> false
    xguide --> "population"
    yguide --> "frequency"
    seriestype --> :histogram
    label --> ids'

    counts
end

function freqplot(result::PartialHistory, tval::AbstractFloat, select)
    t      = result.t
    data   = result.data
    id2ind = result.id2ind

    index   = findfirst(t, tval)
    ix, ids = extract_index_ids(id2ind, select)

    d1 = length(ids)
    d2 = size(data, 3)
    counts = transpose(reshape(data[ix, index, :], d1, d2))

    return counts, ids
end

@userplot SampleTrajectory

@recipe function plot(st::SampleTrajectory; nrlz=1)
    t, x, id = trajplot(st.args..., nrlz)

    # global attributes
    legend -->  true
    grid   --> false
    xguide --> "time"
    yguide --> "population"
    title --> id

    seriestype --> :steppre
    label --> transpose([ "realization $(i)" for i in 1:nrlz ])

    t, x
end

function trajplot(result::PartialHistory, select, n)
    t      = result.t
    data   = result.data
    id2ind = result.id2ind

    ix, ids = extract_index_ids(id2ind, select)

    x = result[ids[1]][:, 1:n]

    return t, x, ids[1]
end

@userplot PhaseTrajectory

@recipe function plot(pt::PhaseTrajectory; nrlz=1)
    x, y, x_id, y_id = phaseplot(pt.args..., nrlz)

    # global attributes
    legend -->  true
    grid   --> false
    xguide --> x_id
    yguide --> y_id

    for i in 1:nrlz
        @series begin
            seriestype --> :path
            label --> "realization $(i)"
            x := x[:, i]
            y := y[:, i]
            ()
        end
    end
end

function phaseplot(result::PartialHistory, x_id, y_id, n)
    t      = result.t
    data   = result.data
    id2ind = result.id2ind

    ix, ids = extract_index_ids(id2ind, [x_id, y_id])

    d1 = length(t)
    d2 = n

    x = result[ids[1]][:, 1:n]
    y = result[ids[2]][:, 1:n]

    return (x, y, ids...)
end
