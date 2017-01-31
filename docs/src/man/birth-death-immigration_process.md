# Birth-Death-Immigration Process

```@setup kendall
using BioSimulator
```

*Kendall's process* is a birth-death-immigration process describing the dynamics of a population using a continuous-time Markov chain. Individuals in the population behave as particles that reproduce at a rate $\alpha$, decay at a rate $\mu$, and immigrate into the population at a rate $\nu$.

```@example kendall
  model = Network("Kendall's Process")

  model <= Species("X", 5)

  model <= Reaction("birth", 2.0, "X --> X + X")
  model <= Reaction("death", 1.0, "X --> 0")
  model <= Reaction("immigration", 0.5, "0 --> X")

  result = simulate(model, algorithm=SSA, time=4.0, epochs=40, trials=1000)

  plot(
    meantrajectory(result),
    freqhistogram(result, 4.0)
  )

  savefig("kendall.svg"); nothing # hide
```
![](kendall.svg)
