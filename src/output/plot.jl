"""
```
plot_results(data::SimulationOutput, ptype::Symbol=:mean; args...)
```

Plot species data according to `ptype`.

### Arguments
- `data`:  a `SimData` object
- `ptype`: The type of graphic to produce. Using `:mean` produces mean trajectories of multiple species, `:sample` plots individual trajectories for a particular species, and `dist` plots a frequency distribution for multiple species at a particular time.
### Optional Arguments
- `ids`: Identifiers for the species to include in a plot. Used by `:mean` and `:dist`.
- `id`: An identifier for a species to include in a plot. Used by `:sample`.
- `nsamples`: The number of trajectories to inclue in a plot. Used by `:sample`.
- `t`: A time used to subset simulation data. Used by `:dist`.
"""
function plot_results(Xt_history::SimData, ptype::Symbol=:mean; args...)
    ptype == :sample ? sampleplot(Xt_history; args...) :
    ptype == :mean   ? meanplot(Xt_history; args...) :
    ptype == :dist   ? distplot(Xt_history; args...) :
    throw(ArgumentError("Unrecognized plot type $(ptype)"))
end

##### dataframe methods #####
function meandf(Xt_history::SimData, ids::Vector{Symbol})
    nrlz = length(Xt_history.rlz)
    ndat = length(Xt_history.rlz[1].t)
    nsps = length(ids)

    dfm  = zeros(Float64, ndat, nsps)
    dfs  = zeros(Float64, ndat, nsps)
    t    = Xt_history.rlz[1].t
    lbl  = Array(Symbol, ndat, nsps)

    j = 1
    for id in ids
        k = Xt_history.id2ind[id]
        Xtk = [ Xt_history.rlz[n].hist[k][i] for n in 1:nrlz, i in 1:ndat ]
        dfm[:,j] = mean(Xtk, 1)[:]
        dfs[:,j] = std(Xtk, 1)[:]
        lbl[:,j] = Symbol[ id for i = 1:ndat]
        j += 1
    end

    df = DataFrame(time = repeat(t, outer=[nsps]),
                    mean = dfm[:],
                    std = dfs[:])
    df[:min] = max(0, df[:mean] - df[:std])
    df[:max] = df[:mean] + df[:std]
    df[:species] = lbl[:]

    return df
end

function sampledf(Xt_history::SimData, nsamples::Int, id::Symbol)
    nrlz = length(Xt_history.rlz)
    ndat = length(Xt_history.rlz[1].t)
    nsps = length(Xt_history.rlz[1].hist)

    #ind = rand(1:nrlz, nsamples)
    t   = Xt_history.rlz[1].t
    k   = get(Xt_history.id2ind, id, 0)
    x   = [ Xt_history.rlz[n].hist[k][i] for i in 1:ndat, n in 1:nsamples ]
    lbl = map(string, 1:nsamples)

    df = DataFrame(time = repeat(t, outer=[nsamples]),
        species_count = x[:],
        iteration = repeat(lbl, inner=[ndat]))
end

function distdf(Xt_history::SimData, ids::Vector{Symbol}, t::Float64)
    nrlz = length(Xt_history.rlz)
    nsps = length(ids)
    i    = findin(Xt_history.rlz[1].t, t)[1]

    Xt  = zeros(Int, nrlz, nsps)
    lbl = map(string, ids)

    j = 1
    for id in ids
        k = Xt_history.id2ind[id]
        Xt[:,j] = Int[ Xt_history.rlz[n].hist[k][i] for n in 1:nrlz ]
        j += 1
    end

    return df = DataFrame(species_count = Xt[:],
        species = repeat(lbl, inner=[nrlz]))
end

##### plotting engines #####
function meanplot(Xt_history::SimData; ids::Vector{Symbol}=Symbol[], na...)
    if isempty(ids)
        error("Must specify species to plot.")
    end

    df = meandf(Xt_history, ids)

    plot(df, x=:time, y=:mean, ymin=:min, ymax=:max, color=:species,
    Geom.line, Geom.point, Geom.errorbar, Theme(major_label_font_size=12pt,minor_label_font_size=12pt)
    )
end

function sampleplot(Xt_history::SimData; nsamples::Integer=1, id::Symbol=:notset, na...)
    if id == :notset
        error("Must specify a species to sample.")
    end

    df = sampledf(Xt_history, nsamples, id)

    plot(df, x=:time, y=:species_count, color=:iteration, Geom.line, Geom.point, Theme(major_label_font_size=12pt,minor_label_font_size=12pt))
end

function distplot(Xt_history::SimData; ids::Vector{Symbol}=Symbol[], t::Float64=-0.0, na...)
    if isempty(ids) || t == -0.0
        error("Must specify species and sample time.")
    end
    df = distdf(Xt_history, ids, t)

    plot(df, x=:species_count, color=:species, Geom.histogram, Theme(major_label_font_size=12pt,minor_label_font_size=12pt))
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
