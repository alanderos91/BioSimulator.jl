immutable MeanTrajectory{T}
    ids :: Vector{T}
end

@recipe function f(mt::MeanTrajectory, h::PartialHistory)
    ids = mt.ids

    x = h.t
    y = h.data
    id2ind = h.id2ind

    legend := true
    grid   := false
    xguide := "time"
    yguide := "mean"

    Xt_mean = transpose(reshape(mean(y, 3), size(y, 1), size(y, 2)))

    for i in eachindex(ids)
        @series begin
            seriestype := :path
            label := ids[i]
            x, Xt_mean[:, id2ind[ids[i]]]
        end
    end
end


@recipe function f(h::PartialHistory)
    # error checking

    x = h.t
    y = h.data
    id2ind = h.id2ind

    legend := true
    grid   := false
    xguide := "time"
    yguide := "mean"
    label  := transpose(map(string, keys(id2ind)))

    Xt_mean = transpose(reshape(mean(y, 3), size(y, 1), size(y, 2)))

    @series begin
        seriestype := :path
        x, Xt_mean
    end
end

immutable Histogram{T,S}
    ids :: Vector{T}
    t   :: S
end

@recipe function f(hg::Histogram, h::PartialHistory)
    ids = hg.ids
    t   = hg.t

    id2ind = h.id2ind

    x = transpose(h[t])

    legend := true
    grid   := false
    xguide := "population"
    yguide := "frequency"

    for i in eachindex(ids)
        @series begin
            seriestype := :histogram
            label := ids[i]
            x[:, id2ind[ids[i]]]
        end
    end
end

immutable SampleTrajectory{T,S}
    ids      :: Vector{T}
    nsamples :: S
end

@recipe function f(st::SampleTrajectory, h::PartialHistory)
    ids = st.ids
    n   = st.nsamples

    x      = h.data
    id2ind = h.id2ind
    t      = h.t

    legend := true
    grid   := false
    xguide := "time"
    yguide := "population"

    layout := grid(length(ids), 1)

    seriestype := :path

    for i in eachindex(ids)
        for j in 1:n
            @series begin
                subplot := i
                title   := string(ids[i])
                label   := "realization $(j)"
                t, vec(x[id2ind[ids[i]], :, j])
            end
        end
    end
end
