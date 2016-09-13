"""
```
simulate(model::Network, algorithm::Algorithm; sampling_interval::AbstractFloat=1.0, nrlz::Integer=1)
```
### Arguments
- `model`: A `Network` to simulate.
- `algorithm`: An algorithm to simulate the model.

### Optional Arguments
- `sampling_interval`: A parameter specifying how frequently output should be recorded. For example, `sampling_interval=0.1` indicates that output is recorded every 0.1 units of time.
- `nrlz`: The number of Monte Carlo realizations to generate.

"""
function simulate(model::Network, algorithm::Algorithm; sampling_interval::AbstractFloat=1.0, nrlz::Integer=1)

  c = n_species(model)
  d = n_reactions(model)

  species   = species_list(model)
  reactions = reaction_list(model)

  X0, id, id2ind = make_species_vector(species)

  if d <= 8
    r = DenseReactionSystem(reactions, id2ind, c, d)
  else
    r = SparseReactionSystem(reactions, id2ind, c, d)
  end

  t  = end_time(algorithm)

  Xt = deepcopy(X0)

  npts = round(Int, t / sampling_interval + 1)

  output = PartialHistory(length(Xt), npts, nrlz, 0.0, t, id2ind)

  simulate(output, Xt, algorithm, X0, r, nrlz)
end

function simulate(output::PartialHistory, Xt::Vector{Int}, algorithm::Algorithm, X0::Vector{Int}, r::AbstractReactionSystem, nrlz::Integer) # nrlz is encoded in PartialHistory; refactor
  a = propensities(r)

  init!(algorithm, Xt, r)
  for i in 1:nrlz
    # setup
    copy!(Xt, X0)
    update_all_propensities!(r, Xt)
    reset!(algorithm, a)
    interval = 1

    while !done(algorithm)
      interval = update!(output, Xt, get_time(algorithm), interval, i)
      step!(algorithm, Xt, r)

      if islossy(a)
        a.intensity = sum(a)
        a.error_bound = zero(eltype(a))
      end
    end

    interval = update!(output, Xt, get_time(algorithm), interval, i)
  end

  return output
end
