@userplot MeanTrajectory

@recipe function f(mt::MeanTrajectory)

    is_one_arg = length(mt.args) == 1 && typeof(mt.args[1]) <: PartialHistory
    is_two_arg = length(mt.args) == 2 && typeof(mt.args[1]) <: PartialHistory && typeof(mt.args[2]) <: Vector

    if is_one_arg || is_two_arg
        result = mt.args[1]

        t      = result.t
        Xt     = result.data
        id2ind = result.id2ind
    else
        error("Mean Trajectory should be given simulation data and, optionally, a vector of identifiers.")
    end

    if is_one_arg
        ids = collect(keys(id2ind)) # plot every species
    elseif is_two_arg
        ids = map(symbol, mt.args[2])
    else
        error("Identifiers should be strings or symbols.")
    end

    Xt_mean = transpose(reshape(mean(Xt, 3), size(Xt, 1, 2)))
    Xt_std  = transpose(reshape(std(Xt, 3), size(Xt, 1, 2)))

    # global attributes
    legend :=  true
    grid   --> false
    xguide --> "time"
    yguide --> "population mean"
    xlims  --> (t[1], t[end])
    ylims  --> (0.0, Inf)
    fillalpha --> 0.3

    # generate each series
    for i in eachindex(ids)
        index = id2ind[ids[i]]

        @series begin
            seriestype --> :path
            label      --> ids[i]

            x      := t
            y      := Xt_mean[:, index]
            ribbon := Xt_std[:, index]

            ()
        end
    end
end

@userplot FrequencyHistogram

@recipe function f(h::FrequencyHistogram)

    is_two_arg = length(h.args) == 2 && typeof(h.args[1]) <: PartialHistory && typeof(h.args[2]) <: AbstractFloat
    is_thr_arg = length(h.args) == 3 && typeof(h.args[1]) <: PartialHistory && typeof(h.args[2]) <: AbstractFloat && typeof(h.args[3]) <: Vector

    if is_two_arg || is_thr_arg
        result = h.args[1]
        t      = h.args[2]
        id2ind = result.id2ind
    else
        error("Frequency Histogram should be given simulation data, a time, and, optionally, a vector of identifiers.")
    end

    if is_two_arg
        ids = collect(keys(id2ind)) # plot every species
    elseif is_thr_arg
        ids = map(symbol, h.args[3])
    else
        error("Identifiers should be strings or symbols.")
    end

    Xt_hist = transpose(result[t])

    legend :=  true
    grid   --> false
    xguide --> "population"
    yguide --> "frequency"

    for i in eachindex(ids)
        index = id2ind[ids[i]]
        @series begin
            seriestype --> :histogram
            label      --> ids[i]

            Xt_hist[:, index]
        end
    end
end

@userplot SampleTrajectory

@recipe function f(st::SampleTrajectory)

    is_two_arg = length(st.args) == 2 && typeof(st.args[1]) <: PartialHistory && typeof(st.args[2]) <: Integer
    is_thr_arg = length(st.args) == 3 && typeof(st.args[1]) <: PartialHistory && typeof(st.args[2]) <: Integer && typeof(st.args[3]) <: Vector

    if is_two_arg || is_thr_arg
        result = st.args[1]
        n      = st.args[2]

        t      = result.t
        Xt     = result.data
        id2ind = result.id2ind
    else
        error("Histogram should be given simulation data, a time, and, optionally, a vector of identifiers.")
    end

    if is_two_arg
        ids = collect(keys(id2ind)) # plot every species
    elseif is_thr_arg
        ids = map(symbol, st.args[3])
    else
        error("Identifiers should be strings or symbols.")
    end

    legend :=  true
    grid   --> false
    xguide --> "time"
    yguide --> "population"
    layout := grid(length(ids), 1)

    for i in eachindex(ids)
        index = id2ind[ids[i]]

        for j in 1:n
            @series begin
                subplot := i
                title  --> string(ids[i])
                label  --> string("realization ", j)

                x := t
                y := vec(Xt[index, :, j])

                ()
            end
        end
    end
end

@userplot PhaseTrajectory

@recipe function f(pt::PhaseTrajectory)
    result = pt.args[1]
    x_id   = pt.args[2]
    y_id   = pt.args[3]
    n      = pt.args[4]

    x_data = result[x_id]
    y_data = result[y_id]

    legend :=  true
    grid   --> false
    xguide --> "population ($x_id)"
    yguide --> "population ($y_id)"

    for j in 1:n
        @series begin
            seriestype --> :path
            label      --> string("realization ", j)

            x := x_data[:, j]
            y := y_data[:, j]

            ()
        end
    end
end
