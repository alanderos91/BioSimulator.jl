struct Configuration{D,T,M,U}
  coord::Vector{NTuple{D,T}}
  tcode::Vector{Int}
  label::Vector{U}
end

function Configuration(lattice::Lattice{D,T,M,U}) where {D,T,M,U}
  coord, tcode = __simple_copy(lattice)
  label = lattice.types

  Configuration{D,T,M,U}(coord, tcode, label)
end

##### plotting recipes

@recipe function f(config::Configuration{D,T,VonNeumann}) where {D,T}
  # unpack information
  coord = config.coord
  tcode = config.tcode
  label = config.label

  # set plot settings
  seriestype := :scatter
  aspect_ratio --> 1
  markerstrokewidth --> 0
  markershape --> :square
  grid --> nothing

  for t in label
    idx = findall(isequal(last(t)), tcode)

    @series begin
      label --> first(t)
      coord[idx]
    end
  end
end
