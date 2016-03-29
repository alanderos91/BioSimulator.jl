# inspired by Mamba.jl

import Gadfly: Geom, Guide, plot

function plot(data::SimulationOutput, ptype::Symbol; field=:species, args...)
    df = getfield(data, field)
    mdata = metadata(data)

    ptype == :sample ? sampleplot(df; args...) :
    ptype == :mean   ? meanplot(df; args...) :
    ptype == :dist   ? distplot(df; args...) :
    throw(ArgumentError("Unrecognized type $(ptype)"))
end

##### plotting engines #####
function meanplot(df; colnames::Vector{Symbol}=names(df)[2:end], na...)
    temp = flatten(df[:, [:time; collect(colnames)]])
    temp = aggregate(temp, [:time, :species], [mean, std])
    names!(temp, [:time, :species, :mean, :std])
    temp[:min] = max(0, temp[:mean] - temp[:std])
    temp[:max] = temp[:mean] + temp[:std]

    plot(temp,
    x=:time, y=:mean, ymin=:min, ymax=:max, color=:species,
    Geom.line, Geom.point, Geom.errorbar, Theme(major_label_font_size=12pt,minor_label_font_size=12pt)
    )
end

function sampleplot(df; nsamples::Integer=1, col::Symbol=names(df)[end], na...)
    temp = df[:, [:time, col]]
    npts  = length(unique(temp[:time]))

    temp = flatten(temp[1:nsamples * npts, :])
    temp[:iteration] = repeat([string(i) for i in 1:nsamples], outer=[1], inner=[npts])

    plot(temp, x=:time, y=:copynumber, color=:iteration, Geom.line, Geom.point, Theme(major_label_font_size=12pt,minor_label_font_size=12pt))
end

function distplot(df; colnames::Vector{Symbol}=names(df)[2:end], t=df[:time][end], na...)
    temp = df[ df[:time] .== t, :]
    temp = flatten(temp[:, [:time; collect(colnames)]])
    plot(temp, x=:copynumber, color=:species, Geom.histogram, Theme(major_label_font_size=12pt,minor_label_font_size=12pt))
end

function flatten(df)
    cols  = names(df)[2:end]
    ncols = length(cols)
    npts  = length(unique(df[:time]))
    itr   = round(Int, length(df[:time]) / npts)

    if ncols == 1
        temp = DataFrame(
        time  = df[:time],
        copynumber = df[cols[1]],
        species = cols[1]
        )
    else
        counts = DataArray(eltype(df[2]), 0)

        for col in cols
            counts = vcat(counts, df[col])
        end

        temp = DataFrame(
        time = repeat(convert(Vector, df[:time]), outer=[ncols]),
        copynumber = counts,
        species = repeat(cols, inner=[npts * itr])
        )
    end

    return temp
end

# write a test for this function; not working
function interpolate(df, final_t, dt, itr)
    t = df[:time]
    npts = round(Int, final_t / dt) + 1
    inds = zeros(Int, itr*npts)
    k = 1
    next_t = 0.0
    for i = 1:length(t)
        while t[i] >= next_t
            inds[k] = i
            next_t = next_t + dt
            k = k + 1
        end

        if next_t > final_t
            next_t = 0.0
        end
    end
    temp = df[inds, :]
    temp[:time] = repeat(collect(linspace(0.0, final_t, npts)), outer=[itr])
    return temp
end
