
````julia
using InteractiveUtils, Plots, BenchmarkTools
using BenchmarkTools, Profile, Random
using NewBioSimulator
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
````





## Old results

Please note that this code is not executed; everything is copied from some old notes.
Importantly, this version did not save simulation data.

#### Direct method with linear search

````julia

@benchmark NewBioSimulator.simulate($state1, $model1, Direct(), 2.0, HasRates)
````


````julia

BenchmarkTools.Trial:
  memory estimate:  209.72 KiB
  allocs estimate:  3411
  --------------
  minimum time:     5.181 ms (0.00% GC)
  median time:      5.531 ms (0.00% GC)
  mean time:        5.656 ms (0.80% GC)
  maximum time:     11.883 ms (46.32% GC)
  --------------
  samples:          884
  evals/sample:     1
````




#### Direct method with bisection

````julia

@benchmark NewBioSimulator.simulate($state1, $model1, Direct(), 2.0, HasSums)
````


````julia

BenchmarkTools.Trial:
  memory estimate:  209.58 KiB
  allocs estimate:  3399
  --------------
  minimum time:     5.177 ms (0.00% GC)
  median time:      5.520 ms (0.00% GC)
  mean time:        5.604 ms (0.80% GC)
  maximum time:     10.697 ms (46.03% GC)
  --------------
  samples:          892
  evals/sample:     1
````




## Test Model


````julia
function example()
  Random.seed!(5357)

  coord = generate_random_2Dpoints((0, 9), (0, 9), 0.8)
  types = rand(["rabbit", "fox"], size(coord, 2))
  state = Lattice(coord, types, nbhood = VonNeumann())

  skeleton = @def_reactions begin
    A + 0 --> 0 + A, α
    B + 0 --> 0 + B, α
    A + 0 --> A + A, β
    B + A --> B + B, γ
    A --> 0, δ1
    B --> 0, δ2
  end α β γ δ1 δ2

  α  = 1.0  # migration rate
  β  = 0.02 # prey reproduction rate
  γ  = 0.02 # predation rate
  δ1 = 0.01 # prey death rate
  δ2 = 0.01 # predator death rate

  param = [α, β, γ, δ1, δ2]

  model = @enumerate_with_sclass skeleton VonNeumann() 2 param

  return state, model
end

state, model = example()
````





## Benchmarks

#### Initialization

````julia
Random.seed!(5357)
@benchmark simulate($state, $model, Direct(), 0.0, HasRates)
````


````
BenchmarkTools.Trial: 
  memory estimate:  153.59 KiB
  allocs estimate:  2014
  --------------
  minimum time:     1.743 ms (0.00% GC)
  median time:      1.767 ms (0.00% GC)
  mean time:        1.831 ms (1.45% GC)
  maximum time:     6.034 ms (66.66% GC)
  --------------
  samples:          2727
  evals/sample:     1
````





#### Simulation
````julia
Random.seed!(5357)
@benchmark simulate($state, $model, Direct(), 2.0, HasRates)
````


````
BenchmarkTools.Trial: 
  memory estimate:  196.13 KiB
  allocs estimate:  3081
  --------------
  minimum time:     2.901 ms (0.00% GC)
  median time:      3.466 ms (0.00% GC)
  mean time:        3.525 ms (1.20% GC)
  maximum time:     9.882 ms (67.89% GC)
  --------------
  samples:          1418
  evals/sample:     1
````





## Profiling

````julia
function run_profile_expirement(state, model)
  simulate(state, model, Direct(), 2.0, HasRates)
  return nothing
end

Profile.init(n = 10^7, delay = 0.001)

# run once to compile, then do the real run
@profile run_profile_expirement(state, model)
Profile.clear()
@profile run_profile_expirement(state, model)

Profile.print(format = :flat, sortedby = :count)
````


````
Count File                        Line Function                           
    
     1 ./deepcopy.jl                 91 _deepcopy_array_t(::Any, ::Type, ::
I...
     1 ./deepcopy.jl                 28 deepcopy                           
    
     1 ./deepcopy.jl                 67 deepcopy_internal(::Any, ::IdDict{A
n...
     1 ./deepcopy.jl                 78 deepcopy_internal                  
    
     1 ./operators.jl               185 !=                                 
    
     1 ./operators.jl                83 ==                                 
    
     1 ...mpleClassEnumeration.jl   294 sample_class_additions!(::NewBioSim
u...
     1 ...mpleClassEnumeration.jl   363 update_classes_neighbors!(::NewBioS
i...
     1 ...mpleClassEnumeration.jl   265 update_classes_particle!(::NewBioSi
m...
     1 .../cCmTD/src/model/ips.jl   213 execute_pairwise!(::NewBioSimulator
....
     1 .../cCmTD/src/model/ips.jl   217 execute_pairwise!(::NewBioSimulator
....
     1 ...r/cCmTD/src/simulate.jl    66 simulate!(::NewBioSimulator.ExactSi
m...
     1 ...r/cCmTD/src/simulate.jl    87 update!                            
    
     1 ...TD/src/state/Lattice.jl    61 __copy(::NewBioSimulator.Lattice{2,
I...
     2 .../cCmTD/src/model/ips.jl   163 execute_jump!(::NewBioSimulator.Lat
t...
     2 .../cCmTD/src/model/ips.jl   149 simulate!(::NewBioSimulator.ExactSi
m...
     3 ./none                         2 run_profile_expirement(::NewBioSimu
l...
     3 ...r/cCmTD/src/simulate.jl    29 simulate                           
    
     4 ./boot.jl                    328 eval                               
    
     4 ./boot.jl                    326 include                            
    
     4 ./client.jl                  436 _start()                           
    
     4 ./client.jl                  267 exec_options(::Base.JLOptions)     
    
     4 ./loading.jl                1038 include_relative(::Module, ::String
)   
     4 ./sysimg.jl                   29 include(::Module, ::String)        
    
     4 ...imulator/notes/build.jl    14 top-level scope                    
    
     4 ...eave/L0NNS/src/Weave.jl   121 #weave#16(::String, ::Symbol, ::Str
i...
     4 .../Weave/L0NNS/src/run.jl    94 #run#29(::String, ::Symbol, ::Strin
g...
     4 .../Weave/L0NNS/src/run.jl   230 capture_output(::Expr, ::Module, ::
B...
     4 .../Weave/L0NNS/src/run.jl   289 eval_chunk(::Weave.CodeChunk, ::Wea
v...
     4 .../Weave/L0NNS/src/run.jl   130 run_chunk(::Weave.CodeChunk, ::Weav
e...
     4 .../Weave/L0NNS/src/run.jl   208 run_code(::Weave.CodeChunk, ::Weave
....
````


