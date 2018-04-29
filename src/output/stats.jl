struct Trajectory
  id      :: String
  t_index :: LinSpace{Float64}
  val     :: Vector{Int}
end

function Trajectory(output :: SimData, species, trial)
  data    = output[species]
  t_index = output.t_index

  val = reshape(data[:, trial], length(t_index))

  return Trajectory(string(species), t_index, val)
end

function Base.show(io :: IO, mt :: Trajectory)
  println(summary(mt))
  println(io, "t", " ", collect(mt.t_index))
  print(io, mt.id, " ", mt.val)
end

struct MeanTrajectory
  id       :: String
  t_index  :: LinSpace{Float64}
  mean_val :: Vector{Float64}
  std_val  :: Vector{Float64}
end

function MeanTrajectory(output :: SimData, species)
  data    = output[species]
  t_index = output.t_index
  n, m    = size(data)

  mean_val = reshape(mean(data, 2), n)
  std_val  = reshape(std(data, 2, mean=mean_val), n)

  return MeanTrajectory(string(species), t_index, mean_val, std_val)
end

function Base.show(io :: IO, mt :: MeanTrajectory)
  println(summary(mt))
  println(io, "t", " ", collect(mt.t_index))
  print(io, mt.id, " ", mt.mean_val)
end

struct Histogram
  id  :: String
  t   :: Float64
  val :: Vector{Int}
end

function Histogram(output :: SimData, species, t)
  data    = output[species]
  t_index = output.t_index

  row = findfirst(t_index, t)
  val = reshape(data[row, :], size(data, 2))

  return Histogram(string(species), t, val)
end

function Base.show(io :: IO, mt :: Histogram)
  println(summary(mt))
  println(io, "t", "=", mt.t)
  print(io, mt.id, " ", mt.val)
end
