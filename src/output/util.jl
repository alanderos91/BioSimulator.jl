import Base.show

immutable SimulationOutput
    species_data::DataFrame
    propensity_data::DataFrame
end

get_species_data(so::SimulationOutput)    = so.species_data
get_propensity_data(so::SimulationOutput) = so.propensity_data

function compile_data(overseer)
    species_data    = DataFrame()
    propensity_data = DataFrame()

    t_observer  = overseer.t_observer
    s_observers = overseer.s_observers
    r_observers = overseer.r_observers

    if !isempty(t_observer.states)
        species_data[t_observer.id]    = t_observer.states
        propensity_data[t_observer.id] = t_observer.states
    end

    for o in s_observers
        species_data[o.id] = o.states
    end

    for o in r_observers
        propensity_data[o.id] = o.states
    end

    return SimulationOutput(species_data, propensity_data)
end

function Base.show(io::IO, x::SimulationOutput)
    @printf io "[species data] %d x %d\n" size(x.species_data, 1) size(x.species_data, 2)
    @printf io "[propensity data] %d x %d\n" size(x.propensity_data, 1) size(x.propensity_data, 2)
end

function plot_species_timeseries(so::SimulationOutput)
    data = so.species_data
    df = aggregate(data, :Time, [mean, std])
    df_summary = DataFrame(Time=Float64[], Mean=Float64[], Min=Float64[], Max=Float64[], Species=UTF8String[])

    for col in names(data)
      if col == :Time; continue; end
      col_mean = symbol(col, "_mean")
      col_std =  symbol(col, "_std")

      temp = DataFrame(
          Time    = df[:Time],
          Mean    = df[col_mean],
          Min     = max(0, df[col_mean] - df[col_std]),
          Max     = df[col_mean] + df[col_std],
          Species = string(col)
      )
      df_summary = vcat(df_summary, temp)
    end
    return plot(df_summary, x=:Time, y=:Mean, ymin=:Min, ymax=:Max, color=:Species, Geom.line, Geom.point, Geom.errorbar)
end
