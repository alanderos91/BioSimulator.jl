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

##### make sure SamplePath works for Configuration
SamplePath(xs::AbstractVector{T}, ts, dims::NTuple{N}) where {T <: Configuration, N} = SamplePath{T, N, typeof(xs), typeof(ts)}(xs, ts)

SamplePath(xs::AbstractVector{T}, ts::AbstractVector) where T <: Configuration = SamplePath(xs, ts, (1, 1))

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

@recipe function f(config::Configuration{2,T,Hexagonal}) where T
  # unpack information
  coord = config.coord
  tcode = config.tcode
  label = config.label

  # set plot settings
  seriestype := :scatter
  aspect_ratio --> 1
  markerstrokewidth --> 0
  markershape --> :hexagon
  grid --> nothing

  A = cos(pi / 3)
  B = sin(pi / 3)

  for t in label
    idx = findall(isequal(last(t)), tcode)
    new_coord = Vector{NTuple{2,Float64}}(undef, length(idx))

    @series begin
      label --> first(t)

      for i in eachindex(idx)
        c = coord[idx[i]]

        x = c[1] + A * c[2]
        y = B * c[2]

        new_coord[i] = (x, y)
      end

      new_coord
    end
  end
end

# display the initial and final configurations
@recipe function f(xw::SamplePath{T}) where T <: Configuration
  cfg_A = xw.u[1]
  cfg_B = xw.u[end]

  layout := 2

  # initial configuration
  @series begin
    subplot := 1
    cfg_A
  end

  # final configuration
  @series begin
    subplot := 2
    cfg_B
  end
end
