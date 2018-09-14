import BioSimulator: stoichiometry, coefficients, scaled_rates, propensities, dependencies, _sort!
m = brusselator(N=10)

c = n_species(m)
d = n_reactions(m)

species   = species_list(m)
reactions = reaction_list(m)

X0, id, id2ind = BioSimulator.make_species_vector(species)

# dense
system = BioSimulator.SparseReactionSystem(reactions, id2ind, c, d)

ix       = sortperm(shuffle!(collect(1:d)))
indexmap = Dict(i => findfirst(isequal(i), ix) for i in eachindex(ix))

V = deepcopy(stoichiometry(system))
U = deepcopy(coefficients(system))
k = deepcopy(scaled_rates(system))
g = deepcopy(dependencies(system))

stoichiometry_map = Dict(indexmap[i] => V[:,i] for i in eachindex(ix))
coefficients_map  = Dict(indexmap[i] => U[:,i] for i in eachindex(ix))
scaled_rates_map  = Dict(indexmap[i] => k[i] for i in eachindex(ix))
dependencies_map  = Dict(indexmap[i] => g[i] for i in eachindex(ix))

_sort!(system, ix)

for i in eachindex(ix)
    @test stoichiometry(system)[:,i] == stoichiometry_map[i]
    @test coefficients(system)[:,i] == coefficients_map[i]
    @test scaled_rates(system)[i] == scaled_rates_map[i]
    @test dependencies(system)[i] == map(x -> indexmap[x], dependencies_map[i])
end
