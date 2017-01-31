# Prokaryotic Auto-Regulatory Gene Network

```@setup gene
using BioSimulator
```

The influence of noise at the cellular level is difficult to capture in deterministic models. Stochastic simulation is appropriate for the study of regulatory mechanisms in genetics, where key species may be present in low numbers.

```@example gene
  function autoreg(;k1=1.0, k1r=10.0, k2=0.01, k3=10.0, k4=1.0, k4r=1.0, k5=0.1, k6=0.01)
    model = Network("auto-regulation")

    model <= Species("gene",   10)
    model <= Species("P2_gene", 0)
    model <= Species("RNA",     0)
    model <= Species("P",       0)
    model <= Species("P2",      0)

    model <= Reaction("repression binding", k1, "gene + P2 --> P2_gene")
    model <= Reaction("reverse repression binding", k1r, "P2_gene --> gene + P2")
    model <= Reaction("transcription", k2, "gene --> gene + RNA")
    model <= Reaction("translation", k3, "RNA --> RNA + P")
    model <= Reaction("dimerization", k4, "P + P --> P2")
    model <= Reaction("dissociation", k4r, "P2 --> P + P")
    model <= Reaction("RNA degradation", k5, "RNA --> 0")
    model <= Reaction("protein degradation", k6, "P --> 0")

    return model
  end

  model = autoreg()

  result = simulate(model, algorithm=SSA, time=1000.0, epochs=500, trials=100)

  plot(
    meantrajectory(result),
    freqhistogram(result, 4.0)
  )

  savefig("gene.svg"); nothing # hide
```
![](gene.svg)
