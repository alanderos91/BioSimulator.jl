using BioSimulator
using DataFrames
using Gadfly

# String Constants
substrate = "S"
enzyme    = "E"
subenz    = "SE"
product   = "P"

r1 = "Binding"
r2 = "Dissociation"
r3 = "Conversion"

# Parameter Values
a1 = 0.00166
a2 = 0.0001
a3 = 0.1

# Initialize network
network = Network("Michaelis-Menten")

# Species Definitions
add_species!(network, substrate, 301)
add_species!(network, enzyme,    120)
add_species!(network, subenz,      0)
add_species!(network, product,     0)

# Reaction Definitions
add_reaction!(network, r1, "a1"); set_parameter!(network, "a1", a1)
add_reactant!(network, r1, substrate)
add_reactant!(network, r1, enzyme)
add_product!( network, r1, subenz)

add_reaction!(network, r2, "a2"); set_parameter!(network, "a2", a2)
add_reactant!(network, r2, subenz)
add_product!( network, r2, substrate)
add_product!( network, r2, enzyme)

add_reaction!(network, r3, "a3"); set_parameter!(network, "a3", a3)
add_reactant!(network, r3, subenz)
add_product!( network, r3, product)
add_product!( network, r3, enzyme)

@time df = simulate(network, tf=50.0, with=:sal, dt=1.0, itr=100);
      df = aggregate(df, :Time, [mean,std])

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

                   Species = repeat([substrate, enzyme, subenz, product], inner=[51]))

plot(new_df, x=:Time, y=:Mean, ymin=:Min, ymax=:Max, color=:Species, Geom.line, Geom.point, Geom.errorbar)
