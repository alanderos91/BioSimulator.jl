using BioSimulator
using DataFrames
using Gadfly

# String Constants
predator = "Predator"
prey = "Prey"

# Parameter Values

# Initialize network
network = Network("Lotka-Volterra")

# Species Definitions
add_species!()

# Reaction Definitions
add_reaction!()

@time df = simulate()
