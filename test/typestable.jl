# using BioSimulator
# import BioSimulator: DenseReactionSystem, compute_propensities!
#
# Xt = [100, 500, 1] # (X, Y, XY)
#
# pre = transpose([
#   0 0 0; # Zero Order
#   1 0 0; # First Order
#   1 1 0; # Second Order, Type A
#   2 0 0  # Second Order, Type B
# ])
#
# post = transpose([
#   1 0 0; # Zero Order
#   0 0 0; # First Order
#   0 0 1; # Second Order, Type A
#   0 1 0  # Second Order, Type B
# ])
#
# v = post - pre
# u = pre
#
# k_id = [:k0, :k1, :k2a, :k2b]
# k    = [0.1, 0.01, 0.25, 0.5]
#
# parameters = Dict{Symbol,Parameter}(
#   k_id[1] => Parameter(k_id[1], k[1]),
#   k_id[2] => Parameter(k_id[2], k[2]),
#   k_id[3] => Parameter(k_id[3], k[3]),
#   k_id[4] => Parameter(k_id[4], k[4])
# )
#
# d = DenseReactionSystem(v, u, k_id)
#
# @code_warntype compute_propensities!(d, Xt, parameters)

using BioSimulator
import BioSimulator: make_species_arr, reaction_system, Model, presimulate!, ODM

Xt, species_id, id2ind = make_species_arr(species_list(Network("")))
rs = reaction_system(reaction_list(Network("")), id2ind)

m = Model(Network("").id, Xt, rs, parameter_list(Network("")), deepcopy(Xt))

@code_warntype presimulate!(ODM(1.0, 100, 1), m)
