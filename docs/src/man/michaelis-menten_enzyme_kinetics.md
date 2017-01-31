# Michaelis-Menten Enzyme Kinetics

```@setup mmek
using BioSimulator
```
Michaelis-Menten enzyme kinetics is a stepwise process combining first- and second order reactions to describe the conversion of a substrate into a product. An enzyme $E$ binds to a substrate $S$ to form a complex $SE$. Conversion does not happen immediately, so $SE$ may revert to its two components or result in a product $P$ and enzyme $E$.

```@example mmek
  model = Network("enzyme kinetics")

  model <= Species("S", 301)
  model <= Species("E", 100)
  model <= Species("SE",  0)
  model <= Species("P",   0)

  model <= Reaction("Binding", 0.00166, "S + E --> SE")
  model <= Reaction("Disociation", 0.0001, "SE --> S + E")
  model <= Reaction("Conversion", 0.1, "SE --> P + E")

  result = simulate(model, algorithm=SSA, time=50.0, epochs=100, trials=1000)

  plot(
    meantrajectory(result),
    freqhistogram(result, 50.0)
  )

  savefig("mmek.svg"); nothing # hide
```
  ![](mmek.svg)
