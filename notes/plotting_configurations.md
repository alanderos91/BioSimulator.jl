
````julia
using InteractiveUtils
using NewBioSimulator
using Plots
````



````
Julia Version 1.1.0
Commit 80516ca202 (2019-01-21 21:24 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin18.2.0)
  CPU: Intel(R) Core(TM) i7-4750HQ CPU @ 2.00GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, haswell)
Environment:
  JULIA_EDITOR = atom  -a
  JULIA_NUM_THREADS = 4
````





## Foxes and rabbits

The following sets up the point and cell type data:

````julia
# example parameters
number_points = 100
number_types = 1
type_list = ["rabbit", "fox"]
xlim = (-10, 10)

# build a list of points
list = Tuple{Int,Int}[]

while length(list) < number_points
  point = tuple(rand(xlim[1]:xlim[2], 2)...)
  if point ∉ list
    push!(list, point)
  end
end

# reinterpret list as a matrix
# this won't be necessary in the future...
coord = zeros(Int, 2, length(list))

for i in eachindex(list)
  coord[1, i] = list[i][1]
  coord[2, i] = list[i][2]
end

# randomly assign a site to be a fox or rabbit
types = [rand(type_list) for i in eachindex(list)]
````





### VonNeumann neighborhoods

````julia
vonlattice = Lattice(coord, types, nbhood = VonNeumann())
````


````
2-D Lattice with VonNeumann neighborhoods
species: 2-element Array{Pair{String,Int64},1}:
 "rabbit" => 2
    "fox" => 3
````



````julia
plot(vonlattice, markersize = 8, grid = nothing)
````


![](figures/plotting_configurations_5_1.svg)



### Hex neighborhoods

````julia
hexlattice = Lattice(coord, types, nbhood = Hexagonal())
````


````
2-D Lattice with Hexagonal neighborhoods
species: 2-element Array{Pair{String,Int64},1}:
 "rabbit" => 2
    "fox" => 3
````



````julia
plot(hexlattice, markersize = 10, grid = nothing)
````


![](figures/plotting_configurations_7_1.svg)
