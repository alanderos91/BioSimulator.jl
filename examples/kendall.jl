using BioSimulator

m = Network("Kendall's Process")

m <= Species("X", 5)

m <= Reaction("Birth", 2.0, :(X --> X + X))
m <= Reaction("Death", 1.0, :(X --> 0))
m <= Reaction("Immigration", 0.5, :(0 --> X))

result = simulate(m, SAL(4.0, 0.125, 100.0, 0.75), sampling_interval=0.1, nrlz=10_000)

plot(result)
