using BioSimulator

# Initialize network
network = Network("Michaelis-Menten")

# Species Definitions
network <= Species("S", 301)
network <= Species("E", 120)
network <= Species("SE",  0)
network <= Species("P",   0)

# Reaction Definitions
network <= Reaction("Binding",     0.00166, :(S + E --> SE))
network <= Reaction("Disociation", 0.0001,  :(SE --> S + E))
network <= Reaction("Conversion",  0.1,     :(SE --> P + E))

# Simulation
srand(5357)
@time result = simulate(network, SAL(50.0, 0.125, 100.0, 0.75), sampling_interval=1.0, nrlz=10_000);

plot(result)
savefig("./docs/mmek_mean.png")

plot(Histogram([:S, :E, :SE, :P], 50.0), result)
savefig("./docs/mmek_hist.png")

plot(SampleTrajectory([:S, :E, :SE, :P], 5), result)
savefig("./docs/mmek_samp.png")

anim = @animate for t in result.t
  plot(Histogram([:P], t), result)
end

gif(anim, "./docs/mmek_anim.gif", fps=5)
