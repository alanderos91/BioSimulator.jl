using BioSimulator

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
network <= Reaction(r1, k1, :(EMT --> MET))
network <= Reaction(r2, k2, :(MET --> EMT))
network <= Reaction(r3, k3, :(MET --> 2*MET))
network <= Reaction(r4, k4, :(MET --> MET + BPP))
network <= Reaction(r5, k5, :(MET --> 2*BPP))
network <= Reaction(r6, k6, :(BPP --> 2*BPP))
network <= Reaction(r7, k7, :(BPP --> 0))
network <= Reaction(r8, k8, :(BPP --> TD))
network <= Reaction(r9, k9, :(TD  --> 0))

@time result = simulate(network, SAL(2*365.0, 0.125, 100.0, 0.75), sampling_interval=1.0, nrlz=10);

plot(result)
