function plot(tr::PopulationTrace, elements::ElementOrFunctionOrLayers...)
  if !isempty(tr)
    name = species(tr)
    tval, pval = toarrays(tr)
    return Gadfly.plot(x=tval, y=pval, elements...)
  else
    # TODO
  end
end
