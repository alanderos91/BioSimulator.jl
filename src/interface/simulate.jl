"""
```
simulate()
```
### Arguments
- ?

### Optional Arguments
- ?

"""
function simulate(model::Network, algorithm::Algorithm; sampling_interval::AbstractFloat=1.0, nrlz::Integer=1)

  c = n_species(model)
  d = n_reactions(model)

  species   = species_list(model)
  reactions = reaction_list(model)

  X0, id, id2ind = make_species_vector(species)

  if d < 8
    r = DenseReactionSystem(reactions, id2ind, c, d)
  else
    r = SparseReactionSystem(reactions, id2ind, c, d)
  end

  simulate(algorithm, X0, r, sampling_interval, nrlz)
end

##### PartialHistory #####
function simulate(algorithm::Algorithm, X0::Vector{Int}, r::AbstractReactionSystem, sampling_interval::AbstractFloat, nrlz::Integer)
  t  = end_time(algorithm)

  Xt = similar(X0)

  npts = round(Int, t / sampling_interval + 1)

  output = PartialHistory(length(Xt), npts, nrlz, 0.0, t)

  simulate(output, Xt, algorithm, X0, r, nrlz)
end

function simulate(output::PartialHistory, Xt::Vector{Int}, algorithm::Algorithm, X0::Vector{Int}, r::AbstractReactionSystem, nrlz::Integer) # nrlz is encoded in PartialHistory; refactor

  init!(algorithm, Xt, r)
  for i in 1:nrlz
    # setup
    copy!(Xt, X0)
    compute_propensities!(r, Xt)
    reset!(algorithm, propensities(r))
    interval = 1

    while !done(algorithm)
      interval = update!(output, Xt, get_time(algorithm), interval, i)
      step!(algorithm, Xt, r)
    end

    interval = update!(output, Xt, get_time(algorithm), interval, i)
  end

  return output
end
