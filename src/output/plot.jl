function plot_mean_timeseries(data, field, colnames...)
    df = getfield(data, field)
    temp = flatten(df[:, [:Time; collect(colnames)]])
    temp = aggregate(temp, [:time, :species], [mean, std])
    names!(temp, [:time, :species, :mean, :std])
    temp[:min] = max(0, temp[:mean] - temp[:std])
    temp[:max] = temp[:mean] + temp[:std]

    plot(temp,
    x=:time, y=:mean, ymin=:min, ymax=:max, color=:species,
    Geom.line, Geom.point, Geom.errorbar
    )
end

function plot_sample_timeseries(data, field, nsamples, col)
    df = getfield(data, field)
    temp = df[:, [:Time; col]]
    npts  = length(unique(df[:Time]))

    temp = flatten(df[1:nsamples * npts, :])
    temp[:iteration] = repeat([string(i) for i in 1:nsamples], outer=[1], inner=[npts])

    plot(temp, x=:time, y=:copy_number, ygroup=:species, color=:iteration,
    Geom.subplot_grid(Geom.line, Geom.point))
end

function plot_histogram(data, field, t, colnames...)
    df = getfield(data, field)
    temp = df[ df[:Time] .== t, :]
    temp = flatten(temp[:, [:Time; collect(colnames)]])
    plot(temp, x=:copy_number, color=:species, Geom.histogram)
end

function flatten(df)
    cols   = names(df)[2:end]
    ncols = length(cols)
    npts  = length(unique(df[:Time]))
    itr    = round(Int, length(df[:Time]) / npts)

    if ncols == 1
        temp = DataFrame(
            time  = df[:Time],
            copynumber = df[cols[1]],
            species = cols[1]
        )
    else
        counts = DataArray(eltype(df[2]), 0)

        for col in cols
            counts = vcat(counts, df[col])
        end

        temp = DataFrame(
            time  = repeat(convert(Vector, df[:Time]), outer=[ncols]),
            copynumber = counts,
            species = repeat(cols, inner=[npts * itr])
        )
    end

    return temp
end
