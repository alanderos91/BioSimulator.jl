# Developers

## The `simulate` Routine

The flexibility of the `simulate` routine allows one to easily implement additional algorithms:

```julia
function serial_simulate(...)

  init!(algorithm, Xt, r)

  for i in 1:max_trials
    copy!(Xt, X0)
    update_all_propensities!(a, r, Xt)
    reset!(algorithm, a)
    interval = 1

    while !done(algorithm)
      interval = update!(output, Xt, get_time(algorithm), interval, i)
      step!(algorithm, Xt, r)
    end

    interval = update!(output, Xt, get_time(algorithm), interval, i)
  end

  return output
end
```

Here is a breakdown of what happens:

- `init!(...)`: This performs any initial computations and data structure setup to avoid allocating memory within a loop.
- `for i in 1:max_trials`: This loop iterates over independent simulations.
- `copy!(...)`: Reset the populations to their initial values.
- `reset!(...)`: Reset algorithm variables and data structures.
- `while !done(...)`: The main loop in the algorithm.
- `update!(...)`: Sample the state vector and record it as needed.
- `step!(...)`: Execute one step in the algorithm.

Thus, one simply needs to implement the `init!`, `reset!`, and `step!` functions for an algorithm.

## The `Algorithm` Interface

```@docs
  Algorithm
```

### The `ExactMethod` Interface

```@docs
  ExactMethod
```

### The `TauLeapMethod` Interface

```@docs
  TauLeapMethod
```
