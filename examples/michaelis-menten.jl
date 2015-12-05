using BioSimulator
using DataFrames
using Gadfly

# Initialize network
network = Network("Michaelis-Menten")

# Species Definitions
network <= Species(:S, 301)
network <= Species(:E, 120)
network <= Species(:SE,  0)
network <= Species(:P,   0)

# Reaction Definitions
network <= Reaction(:Binding,     :a1, r=(:S=>1, :E=>1), p=(:SE=>1))
network <= Reaction(:Disociation, :a2, r=(:SE=>1),       p=(:S=>1, :E=>1))
network <= Reaction(:Conversion,  :a3, r=(:SE=>1),       p=(:P=>1, :E=>1))

# Parameter Definitions
network <= parameter(:a1, 0.00166)
network <= parameter(:a2, 0.0001)
network <= parameter(:a3, 0.1)

# Simulation
srand(5357)
@time df = simulate(network, tf=50.0, with=:sal, dt=1.0, itr=100);
      df = aggregate(df, :Time, [mean,std])

# Convert default dataframe
new_df = DataFrame(Time    = repeat(convert(Vector, df[:Time]), outer=[4]),
                   Mean    = vcat(df[:S_mean], df[:E_mean], df[:SE_mean], df[:P_mean]),
                   Std     = vcat(df[:S_std],  df[:E_std],  df[:SE_std],  df[:P_std]),

                   Min     = vcat(df[:S_mean]  - df[:S_std],
                                 df[:E_mean]  - df[:E_std],
                                 df[:SE_mean] - df[:SE_std],
                                 df[:P_mean]  - df[:P_std]),

                   Max     = vcat(df[:S_mean] + df[:S_std],
                                 df[:E_mean]  + df[:E_std],
                                 df[:SE_mean] + df[:SE_std],
                                 df[:P_mean]  + df[:P_std]),

                   Species = repeat([:S, :E, :SE, :P], inner=[51]))

# Plot results
plot(new_df, x=:Time, y=:Mean, ymin=:Min, ymax=:Max, color=:Species, Geom.line, Geom.point, Geom.errorbar)
