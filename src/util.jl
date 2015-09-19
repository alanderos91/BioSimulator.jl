function set_metadata!(md, alg, t_final, itr, output, stepsize)
  md["algorithm"] = alg
  md["iterations"] = itr
  md["duration"] = t_final
  if output == :explicit
    md["output"] = "Explicit"
  elseif output == :fixed && stepsize > 0.0
    md["output"] = "Fixed-Interval"
    md["stepsize"] = stepsize
  end
end
