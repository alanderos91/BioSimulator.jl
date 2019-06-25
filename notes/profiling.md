
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
  memory estimate:  99.73 KiB
  allocs estimate:  971
  --------------
  minimum time:     749.701 μs (0.00% GC)
  median time:      906.217 μs (0.00% GC)
  mean time:        907.704 μs (2.75% GC)
  maximum time:     72.333 ms (98.71% GC)
  --------------
  samples:          5497
  evals/sample:     1
````





#### Simulation
````julia
Random.seed!(5357)
@benchmark simulate($state, $model, Direct(), 2.0, HasRates)
````


````
BenchmarkTools.Trial: 
  memory estimate:  1.87 MiB
  allocs estimate:  5561
  --------------
  minimum time:     2.264 ms (0.00% GC)
  median time:      2.968 ms (0.00% GC)
  mean time:        3.175 ms (7.59% GC)
  maximum time:     71.224 ms (95.22% GC)
  --------------
  samples:          1573
  evals/sample:     1
````





## Profiling

````julia
function run_profile_expirement(state, model)
  simulate(state, model, Direct(), 2.0, HasRates)
  return nothing
end

Profile.init(n = 10^7, delay = 1e-5)

# run once to compile, then do the real run
@profile run_profile_expirement(state, model)
Profile.clear()
@profile run_profile_expirement(state, model)

Profile.print(format = :flat, sortedby = :count)
````


````
Count File                        Line Function                           
    
     1 ./abstractarray.jl           434 checkbounds                        
    
     1 ./abstractarray.jl           449 checkbounds                        
    
     1 ./abstractarray.jl           506 checkindex                         
    
     1 ./abstractarray.jl           519 checkindex                         
    
     1 ./abstractarray.jl          1711 cmp(::Array{Int64,1}, ::Array{Int64
,1})
     1 ./abstractarray.jl           635 empty                              
    
     1 ./abstractarray.jl           573 similar                            
    
     1 ./array.jl                   729 _zip_iterate_some                  
    
     1 ./array.jl                   274 copyto!(::Array{Int64,1}, ::Int64, 
:...
     1 ./array.jl                   286 copyto!                            
    
     1 ./array.jl                   668 grow_to!(::Array{Int64,1}, ::Base.G
e...
     1 ./array.jl                   155 istracked                          
    
     1 ./array.jl                   729 iterate                            
    
     1 ./array.jl                   199 length                             
    
     1 ./array.jl                   729 mt_pop!                            
    
     1 ./array.jl                   317 similar                            
    
     1 ./array.jl                   812 update!                            
    
     1 ./boot.jl                    421 Type                               
    
     1 ./broadcast.jl               544 _broadcast_getindex                
    
     1 ./broadcast.jl               550 _broadcast_getindex                
    
     1 ./broadcast.jl               575 _getindex                          
    
     1 ./deepcopy.jl                 85 _deepcopy_array_t(::Any, ::Type, ::
I...
     1 ./deepcopy.jl                 88 _deepcopy_array_t(::Any, ::Type, ::
I...
     1 ./deepcopy.jl                 56 deepcopy_internal(::Any, ::IdDict{A
n...
     1 ./deepcopy.jl                 75 deepcopy_internal(::Array{Int64,1},
 ...
     1 ./float.jl                   399 *                                  
    
     1 ./float.jl                   565 hash                               
    
     1 ./hashing.jl                  18 hash                               
    
     1 ./int.jl                      52 -                                  
    
     1 ./int.jl                     444 >>>                                
    
     1 ./int.jl                     452 >>>                                
    
     1 ./int.jl                     136 abs                                
    
     1 ./int.jl                      96 flipsign                           
    
     1 ./int.jl                     444 searchsortedlast                   
    
     1 ./iterators.jl               350 _zip_iterate_some                  
    
     1 ./iterators.jl               352 _zip_iterate_some                  
    
     1 ./iterators.jl               427 iterate                            
    
     1 ./iterators.jl               429 iterate                            
    
     1 ./iterators.jl               169 pairs                              
    
     1 ./iterators.jl               224 pairs                              
    
     1 ./multidimensional.jl        642 __simple_copy(::NewBioSimulator.Lat
t...
     1 ./multidimensional.jl        641 _getindex                          
    
     1 ./multidimensional.jl        654 _unsafe_getindex(::IndexLinear, ::A
r...
     1 ./operators.jl               949 in                                 
    
     1 ./ordering.jl                 51 lt                                 
    
     1 ./pointer.jl                  65 unsafe_copyto!                     
    
     1 ./sort.jl                    293 #searchsortedfirst#4               
    
     1 ./sort.jl                    293 #searchsortedfirst#4               
    
     1 ./sort.jl                    293 #searchsortedfirst#4(::Function, ::
F...
     1 ./sort.jl                    176 sample_class_removals!(::NewBioSimu
l...
     1 ./sort.jl                    192 searchsorted(::Array{NewBioSimulato
r...
     1 ./sort.jl                    177 searchsortedfirst                  
    
     1 ./sort.jl                    191 searchsortedlast                   
    
     1 ./sort.jl                    540 sort!(::Array{NewBioSimulator.Site{
2...
     1 ./sort.jl                    546 sort!(::Array{NewBioSimulator.Site{
2...
     1 ./sort.jl                    556 sort!(::Array{NewBioSimulator.Site{
2...
     1 ./sort.jl                    176 update_classes_neighbor_change!(::N
e...
     1 ./sort.jl                    181 update_classes_neighbor_change!(::N
e...
     1 ./stream.jl                  394 close(::Base.PipeEndpoint)         
    
     1 ./sysimg.jl                   18 get_ptype                          
    
     1 ./sysimg.jl                   19 setproperty!                       
    
     1 ./tuple.jl                   339 <                                  
    
     1 ./tuple.jl                   366 isless                             
    
     1 .../abstract_algorithms.jl    34 generate_jump                      
    
     1 .../abstract_algorithms.jl    23 update_jump_rates!                 
    
     1 ...rc/algorithms/direct.jl    10 update!                            
    
     1 ...algorithms/ssa_utils.jl    10 update_jump_rates!                 
    
     1 ...mpleClassEnumeration.jl   182 derive_nbclass(::NewBioSimulator.Sa
m...
     1 ...mpleClassEnumeration.jl   158 get_nbclass_index                  
    
     1 ...mpleClassEnumeration.jl   200 get_new_nbclass(::NewBioSimulator.S
a...
     1 ...mpleClassEnumeration.jl   277 sample_class_removals!(::NewBioSimu
l...
     1 ...mpleClassEnumeration.jl   282 sample_class_removals!(::NewBioSimu
l...
     1 ...mpleClassEnumeration.jl   230 sample_from_class                  
    
     1 ...mpleClassEnumeration.jl   427 sample_neighbor(::NewBioSimulator.S
i...
     1 ...mpleClassEnumeration.jl   452 sample_neighbor(::Base.Generator{Ar
r...
     1 ...mpleClassEnumeration.jl   317 update_classes_neighbor_change!(::N
e...
     1 ...mpleClassEnumeration.jl   372 update_classes_neighbors!(::NewBioS
i...
     1 ...mpleClassEnumeration.jl   375 update_classes_neighbors!(::NewBioS
i...
     1 ...mpleClassEnumeration.jl   388 update_classes_neighbors!(::NewBioS
i...
     1 ...mpleClassEnumeration.jl   393 update_classes_neighbors!(::NewBioS
i...
     1 ...mpleClassEnumeration.jl   397 update_classes_neighbors!(::NewBioS
i...
     1 ...mpleClassEnumeration.jl   400 update_classes_neighbors!(::NewBioS
i...
     1 ...mpleClassEnumeration.jl   411 update_classes_neighbors!(::NewBioS
i...
     1 .../FbL9C/src/model/ips.jl   177 execute_jump!(::NewBioSimulator.Lat
t...
     1 .../FbL9C/src/model/ips.jl   189 execute_pairwise!(::NewBioSimulator
....
     1 .../FbL9C/src/model/ips.jl   205 execute_pairwise!(::NewBioSimulator
....
     1 .../FbL9C/src/model/ips.jl   210 execute_pairwise!(::NewBioSimulator
....
     1 .../FbL9C/src/model/ips.jl   256 rate                               
    
     1 .../FbL9C/src/model/ips.jl   270 rate                               
    
     1 ...output/Configuration.jl     8 Type                               
    
     1 ...r/FbL9C/src/simulate.jl    60 build_output                       
    
     1 ...r/FbL9C/src/simulate.jl    95 initialize_datastructs!(::NewBioSim
u...
     1 ...r/FbL9C/src/simulate.jl    97 initialize_datastructs!(::NewBioSim
u...
     1 ...r/FbL9C/src/simulate.jl   108 initialize_datastructs!(::NewBioSim
u...
     1 ...r/FbL9C/src/simulate.jl   112 initialize_datastructs!(::NewBioSim
u...
     1 ...r/FbL9C/src/simulate.jl   119 initialize_datastructs!(::NewBioSim
u...
     1 ...r/FbL9C/src/simulate.jl    24 simulate                           
    
     1 ...src/simulators/exact.jl    56 generate_next_jump!                
    
     1 ...src/simulators/exact.jl    45 step!                              
    
     1 ...src/simulators/exact.jl    48 step!                              
    
     1 ...9C/src/state/Lattice.jl   129 add_neighbor!(::NewBioSimulator.Lat
t...
     1 ...9C/src/state/Lattice.jl   170 add_neighbor!(::NewBioSimulator.Lat
t...
     1 ...9C/src/state/Lattice.jl   122 add_site!(::NewBioSimulator.Lattice
{...
     1 ...9C/src/state/Lattice.jl   123 add_site!(::NewBioSimulator.Lattice
{...
     1 ...9C/src/state/Lattice.jl   264 build_neighborhoods!(::NewBioSimula
t...
     1 ...9C/src/state/Lattice.jl   267 build_neighborhoods!(::NewBioSimula
t...
     1 ...FbL9C/src/state/site.jl    28 Type                               
    
     1 ...FbL9C/src/state/site.jl    32 Type                               
    
     1 ...FbL9C/src/state/site.jl    97 build_composition!(::Array{Int64,1}
,...
     1 ...FbL9C/src/state/site.jl    15 change_neighbor_class!             
    
     1 ...FbL9C/src/state/site.jl    51 change_neighbor_class!             
    
     1 ...FbL9C/src/state/site.jl   149 sort_by_x_2D!                      
    
     1 ...FbL9C/src/state/site.jl   154 sort_by_y_2D!                      
    
     1 ...FbL9C/src/state/site.jl    50 transform!                         
    
     1 ...1Rb/src/dict_support.jl     6 hashindex                          
    
     1 ...1Rb/src/ordered_dict.jl   353 getindex                           
    
     1 ...1Rb/src/ordered_dict.jl   228 ht_keyindex(::OrderedCollections.Or
d...
     1 ...1Rb/src/ordered_dict.jl   231 ht_keyindex(::OrderedCollections.Or
d...
     1 ...1Rb/src/ordered_dict.jl   242 ht_keyindex(::OrderedCollections.Or
d...
     1 .../Weave/L0NNS/src/run.jl   246 capture_output(::Expr, ::Module, ::
B...
     1 ...v1.1/Random/src/RNGs.jl   315 rand                               
    
     1 ...v1.1/Random/src/RNGs.jl   306 rand_inbounds                      
    
     1 ...v1.1/Random/src/RNGs.jl   310 rand_inbounds                      
    
     1 ...andom/src/generation.jl   119 rand                               
    
     1 ....1/Random/src/normal.jl    97 randexp                            
    
     1 ....1/Random/src/normal.jl    98 randexp(::Random.MersenneTwister)  
    
     2 ./abstractarray.jl          1706 cmp(::Array{Int64,1}, ::Array{Int64
,1})
     2 ./abstractdict.jl            588 setindex!                          
    
     2 ./array.jl                   669 grow_to!(::Array{Int64,1}, ::Base.G
e...
     2 ./array.jl                   698 grow_to!(::Array{Int64,1}, ::Base.G
e...
     2 ./deepcopy.jl                 63 deepcopy_internal(::Any, ::IdDict{A
n...
     2 ./generator.jl                44 iterate                            
    
     2 ./int.jl                     428 <=                                 
    
     2 ./iterators.jl               342 _zip_iterate_all                   
    
     2 ./iterators.jl               429 grow_to!(::Array{Int64,1}, ::Base.G
e...
     2 ./iterators.jl               332 iterate                            
    
     2 ./operators.jl               185 !=                                 
    
     2 ./operators.jl                83 ==                                 
    
     2 ./sort.jl                    293 #searchsortedfirst#4               
    
     2 ./sort.jl                    694 #sort!#7                           
    
     2 ./sort.jl                    210 searchsorted(::Array{NewBioSimulato
r...
     2 ./sort.jl                    216 searchsorted(::Array{NewBioSimulato
r...
     2 ./sort.jl                    217 searchsorted(::Array{NewBioSimulato
r...
     2 ./sort.jl                    176 searchsortedfirst                  
    
     2 ./sort.jl                    539 sort!                              
    
     2 ./sort.jl                    634 sort!                              
    
     2 ./tuple.jl                   319 hash                               
    
     2 ...mpleClassEnumeration.jl   217 add_class_member!(::Array{Int64,1},
 ...
     2 ...mpleClassEnumeration.jl   285 sample_class_removals!(::NewBioSimu
l...
     2 ...mpleClassEnumeration.jl   337 update_classes_neighbor_change!(::N
e...
     2 ...mpleClassEnumeration.jl   347 update_classes_neighbor_change!(::N
e...
     2 ...mpleClassEnumeration.jl   367 update_classes_neighbors!(::NewBioS
i...
     2 ...mpleClassEnumeration.jl   403 update_classes_neighbors!(::NewBioS
i...
     2 ...r/FbL9C/src/simulate.jl   103 initialize_datastructs!(::NewBioSim
u...
     2 ...r/FbL9C/src/simulate.jl    39 simulate!(::NewBioSimulator.ExactSi
m...
     2 ...FbL9C/src/state/site.jl     2 Type                               
    
     2 ...FbL9C/src/state/site.jl     8 Type                               
    
     2 ...FbL9C/src/state/site.jl    11 get_ptype                          
    
     2 ...FbL9C/src/state/site.jl    37 state                              
    
     2 ....1/Random/src/Random.jl   219 rand                               
    
     3 ./array.jl                   823 _deleteat!                         
    
     3 ./array.jl                  1175 deleteat!                          
    
     3 ./array.jl                   855 push!                              
    
     3 ./boot.jl                    402 findall                            
    
     3 ./int.jl                      53 +                                  
    
     3 ./sort.jl                    293 searchsortedfirst                  
    
     3 ./sort.jl                    545 sort!(::Array{NewBioSimulator.Site{
2...
     3 ...mpleClassEnumeration.jl   216 add_class_member!(::Array{Int64,1},
 ...
     3 ...mpleClassEnumeration.jl   224 rmv_class_member!                  
    
     3 ...mpleClassEnumeration.jl   384 update_classes_neighbors!(::NewBioS
i...
     3 ...9C/src/state/Lattice.jl   109 spawn_new_site                     
    
     3 ...FbL9C/src/state/site.jl    45 get_ptype                          
    
     3 ...1Rb/src/ordered_dict.jl   358 get                                
    
     4 ./abstractarray.jl           617 similar                            
    
     4 ./abstractarray.jl           618 similar                            
    
     4 ./array.jl                   814 _growat!                           
    
     4 ./array.jl                  1152 insert!                            
    
     4 ./boot.jl                    419 Type                               
    
     4 ./boot.jl                    402 materialize                        
    
     4 ./broadcast.jl               196 similar                            
    
     4 ./cartesian.jl                64 macro expansion                    
    
     4 ./deepcopy.jl                 58 deepcopy_internal(::Any, ::IdDict{A
n...
     4 ./multidimensional.jl        642 _getindex                          
    
     4 ./multidimensional.jl        656 _unsafe_getindex(::IndexLinear, ::A
r...
     4 ./multidimensional.jl        662 _unsafe_getindex!                  
    
     4 ./multidimensional.jl        666 macro expansion                    
    
     4 ./multidimensional.jl        671 macro expansion                    
    
     4 ./sort.jl                    291 searchsortedfirst                  
    
     4 ...mpleClassEnumeration.jl   302 sample_class_additions!(::NewBioSim
u...
     4 ...mpleClassEnumeration.jl   269 update_classes_particle!(::NewBioSi
m...
     5 ./abstractarray.jl           927 getindex                           
    
     5 ./abstractdict.jl            594 get                                
    
     5 ./abstractdict.jl             17 haskey                             
    
     5 ./abstractdict.jl            665 in                                 
    
     5 ./boot.jl                    411 Type                               
    
     5 ./sort.jl                    209 searchsorted(::Array{NewBioSimulato
r...
     5 ./tuple.jl                   336 <                                  
    
     5 ...9C/src/state/Lattice.jl    69 __simple_copy(::NewBioSimulator.Lat
t...
     6 ./abstractarray.jl          1714 isless                             
    
     6 ./array.jl                   705 iterate                            
    
     6 ./boot.jl                    402 Type                               
    
     6 ./ordering.jl                 49 lt                                 
    
     6 ./ordering.jl                 50 lt                                 
    
     6 ./sort.jl                    291 searchsortedfirst(::Array{Array{Int
6...
     7 ./broadcast.jl               551 _broadcast_getindex                
    
     7 ./broadcast.jl               578 _broadcast_getindex_evalf          
    
     7 ...r/FbL9C/src/simulate.jl    26 simulate                           
    
     7 ...FbL9C/src/state/site.jl    38 coordinates                        
    
     8 ./broadcast.jl               511 getindex                           
    
     8 ...mpleClassEnumeration.jl   376 update_classes_neighbors!(::NewBioS
i...
     9 ./sort.jl                    178 searchsortedfirst                  
    
     9 .../FbL9C/src/model/ips.jl   213 execute_pairwise!(::NewBioSimulator
....
    10 ./array.jl                   767 setindex!                          
    
    10 ./ordering.jl                 69 #1                                 
    
    10 ./ordering.jl                 52 lt                                 
    
    10 ...mpleClassEnumeration.jl   364 update_classes_neighbors!(::NewBioS
i...
    10 ...mpleClassEnumeration.jl   244 update_classes_particle!(::NewBioSi
m...
    10 .../FbL9C/src/model/ips.jl   214 execute_pairwise!(::NewBioSimulator
....
    10 ...9C/src/state/Lattice.jl    65 __simple_copy(::NewBioSimulator.Lat
t...
    10 ...9C/src/state/Lattice.jl    82 get_site                           
    
    11 ./ordering.jl                 75 ord(::Function, ::Function, ::Bool,
 ...
    11 ./sort.jl                    293 #searchsortedfirst#4(::Function, ::
F...
    11 ./sysimg.jl                   18 getproperty                        
    
    12 ./array.jl                   729 getindex                           
    
    12 ...9C/src/state/Lattice.jl    64 __simple_copy(::NewBioSimulator.Lat
t...
    13 ./sort.jl                    211 searchsorted(::Array{NewBioSimulato
r...
    14 ./broadcast.jl               797 copyto!                            
    
    14 ./broadcast.jl               842 copyto!                            
    
    14 ./broadcast.jl               843 macro expansion                    
    
    14 ./simdloop.jl                 73 macro expansion                    
    
    14 ./sort.jl                    293 #searchsortedfirst#4(::Function, ::
F...
    16 ./array.jl                   812 _growend!                          
    
    16 ./array.jl                   854 push!                              
    
    16 ...mpleClassEnumeration.jl   363 update_classes_neighbors!(::NewBioS
i...
    17 ./array.jl                   693 grow_to!(::Array{Int64,1}, ::Base.G
e...
    17 ...9C/src/state/Lattice.jl    75 istracked                          
    
    18 ./broadcast.jl               773 copy                               
    
    18 ./broadcast.jl               753 materialize                        
    
    19 ...mpleClassEnumeration.jl   370 update_classes_neighbors!(::NewBioS
i...
    21 ./array.jl                   670 grow_to!(::Array{Int64,1}, ::Base.G
e...
    24 ./array.jl                   604 collect                            
    
    25 ./array.jl                  1966 findall                            
    
    26 ...mpleClassEnumeration.jl   206 get_new_nbclass(::NewBioSimulator.S
a...
    27 ./sort.jl                    293 #searchsorted#6                    
    
    27 ./sort.jl                    291 searchsorted                       
    
    28 ...9C/src/state/Lattice.jl    67 __simple_copy(::NewBioSimulator.Lat
t...
    30 .../FbL9C/src/model/ips.jl   218 execute_pairwise!(::NewBioSimulator
....
    38 .../FbL9C/src/model/ips.jl   217 execute_pairwise!(::NewBioSimulator
....
    41 ./deepcopy.jl                 91 _deepcopy_array_t(::Any, ::Type, ::
I...
    43 ./deepcopy.jl                 28 deepcopy                           
    
    43 ./deepcopy.jl                 78 deepcopy_internal(::Array{NewBioSim
u...
    43 ...r/FbL9C/src/simulate.jl    18 simulate                           
    
    43 ...9C/src/state/Lattice.jl    58 copy                               
    
    55 ...output/Configuration.jl     8 update!                            
    
    56 ...r/FbL9C/src/simulate.jl    40 simulate!(::NewBioSimulator.ExactSi
m...
    84 ./deepcopy.jl                 67 deepcopy_internal(::Any, ::IdDict{A
n...
    91 .../FbL9C/src/model/ips.jl   164 execute_jump!(::NewBioSimulator.Lat
t...
    92 .../FbL9C/src/model/ips.jl   150 simulate!(::NewBioSimulator.ExactSi
m...
   150 ...r/FbL9C/src/simulate.jl    29 simulate                           
    
   201 ./boot.jl                    328 eval                               
    
   201 ./none                         2 run_profile_expirement(::NewBioSimu
l...
   201 .../Weave/L0NNS/src/run.jl   230 capture_output(::Expr, ::Module, ::
B...
   202 ./boot.jl                    326 include                            
    
   202 ./client.jl                  436 _start()                           
    
   202 ./client.jl                  267 exec_options(::Base.JLOptions)     
    
   202 ./loading.jl                1038 include_relative(::Module, ::String
)   
   202 ./sysimg.jl                   29 include(::Module, ::String)        
    
   202 ...imulator/notes/build.jl    14 top-level scope                    
    
   202 ...eave/L0NNS/src/Weave.jl   121 #weave#16(::String, ::Symbol, ::Str
i...
   202 .../Weave/L0NNS/src/run.jl    94 #run#29(::String, ::Symbol, ::Strin
g...
   202 .../Weave/L0NNS/src/run.jl   289 eval_chunk(::Weave.CodeChunk, ::Wea
v...
   202 .../Weave/L0NNS/src/run.jl   130 run_chunk(::Weave.CodeChunk, ::Weav
e...
   202 .../Weave/L0NNS/src/run.jl   208 run_code(::Weave.CodeChunk, ::Weave
....
````


