using BioSimulator
using DataFrames
using Gadfly

# Reaction Names
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
network <= Species(:EMT,   200)
network <= Species(:MET,   800)
network <= Species(:BPP, 99000)
network <= Species(:TD,      0)

# Reaction Definitions
network <= Reaction(r1, :k1, r=(:EMT => 1), p=(:MET => 1))
network <= Reaction(r2, :k2, r=(:MET => 1), p=(:EMT => 1))
network <= Reaction(r3, :k3, r=(:MET => 1), p=(:MET => 2))
network <= Reaction(r4, :k4, r=(:MET => 1), p=(:MET => 1, :BPP => 1))
network <= Reaction(r5, :k5, r=(:MET => 1), p=(:BPP => 2))
network <= Reaction(r6, :k6, r=(:BPP => 1), p=(:BPP => 2))
network <= Reaction(r7, :k7, r=(:BPP => 1))
network <= Reaction(r8, :k8, r=(:BPP => 1), p=(:TD => 1))
network <= Reaction(r9, :k9, r=(:TD => 1))

# Parameter Definitions
network <= parameter(:k1, k1)
network <= parameter(:k2, k2)
network <= parameter(:k3, k3)
network <= parameter(:k4, k4)
network <= parameter(:k5, k5)
network <= parameter(:k6, k6)
network <= parameter(:k7, k7)
network <= parameter(:k8, k8)
network <= parameter(:k9, k9)

@time df = simulate(network, tf=10 * 365.0, with=:ssa, dt=1.0, output=Uniform(), itr=100);
df = aggregate(df, :Time, mean);

p1 = plot(df, x=:Time, y=:EMT_mean, Geom.line)
p2 = plot(df, x=:Time, y=:MET_mean, Geom.line)
p3 = plot(df, x=:Time, y=:BPP_mean, Geom.line)
p4 = plot(df, x=:Time, y=:TD_mean,  Geom.line)
;
