using BioSimulator
using DataFrames
using Gadfly

# String Constants
emt = "EMT"
met = "MET"
bpp = "BPP"
td  = "TD"

r1 = "Mesenchymal Transition"
r2 = "Endometrial Transition"
r3 = "Symmetric Self-Renewal (BCSC)"
r4 = "Asymmetric Self-Renewal (BCSC)"
r5 = "Symmetric Differentiation (BCSC)"
r6 = "Symmetric Self-Renewal (BPP)"
r7 = "Death of Progenitor Cells"
r8 = "Differentiation of Progenitor Cells"
r9 = "Death of Terminally Differentiated Cells"

# Parameter Values
k1 = 0.08
k2 = 0.02
k3 = 0.0034
k4 = 0.027
k5 = 0.0034
k6 = 0.31
k7 = 0.015
k8 = 0.3
k9 = 0.01

# Initialize network
network = Network("BCSC Carcinogenesis")

# Species Definitions
add_species!(network, emt, 200)
add_species!(network, met, 800)
add_species!(network, bpp, 99000)
add_species!(network, td,  0)

# Reaction Definitions
add_reaction!(network, r1, "k1"); set_parameter!(network, "k1", k1)
add_reactant!(network, r1, emt)
add_product!( network, r1, met)

add_reaction!(network, r2, "k2"); set_parameter!(network, "k2", k2)
add_reactant!(network, r2, met)
add_product!( network, r2, emt)

add_reaction!(network, r3, "k3"); set_parameter!(network, "k3", k3)
add_reactant!(network, r3, met)
add_product!( network, r3, met, 2)

add_reaction!(network, r4, "k4"); set_parameter!(network, "k4", k4)
add_reactant!(network, r4, met)
add_product!( network, r4, met)
add_product!( network, r4, bpp)

add_reaction!(network, r5, "k5"); set_parameter!(network, "k5", k5)
add_reactant!(network, r5, met)
add_product!( network, r5, bpp, 2)

add_reaction!(network, r6, "k6"); set_parameter!(network, "k6", k6)
add_reactant!(network, r6, bpp)
add_product!( network, r6, bpp, 2)

add_reaction!(network, r7, "k7"); set_parameter!(network, "k7", k7)
add_reactant!(network, r7, bpp)

add_reaction!(network, r8, "k8"); set_parameter!(network, "k8", k8)
add_reactant!(network, r8, bpp)
add_product!( network, r8, td)

add_reaction!(network, r9, "k9"); set_parameter!(network, "k9", k9)
add_reactant!(network, r9, td)

@time df = simulate(network, tf=100.0, with=:sal, dt=1.0, output=Uniform(), itr=10); df = aggregate(df, :Time, mean);

plot(df, x=:Time, y=:EMT_mean, Geom.line)
plot(df, x=:Time, y=:MET_mean, Geom.line)
plot(df, x=:Time, y=:BPP_mean, Geom.line)
plot(df, x=:Time, y=:TD_mean,  Geom.line)
