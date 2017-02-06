var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#BioSimulator.jl-1",
    "page": "Home",
    "title": "BioSimulator.jl",
    "category": "section",
    "text": "A stochastic simulation framework for Julia."
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "Many complex systems in biology are analytically intractable, and dynamical predictions based on deterministic models can be grossly misleading. Stochastic simulation algorithms based on continuous-time Markov chains allow researchers to generate accurate time-evolution trajectories, test the sensitivity of models to key parameters, and quantify frequencies of rare events [6, 11, 15].Situations where stochastic simulation is especially helpful involve:rare events such as extinction and mutation,\nkey molecules present in small numbers,\nrare reactions with dramatic influence on the dynamics of the system\npopulation cycles arising from demographic stochasticity.Examples of such systems include gene expression networks, tumor suppressor pathways, and demographic and ecological systems.BioSimulator.jl aims to provide researchers interested in such phenomena with a fast, reliable, user-friendly, and open-source modeling tool in Julia.This package supports Julia 0.5."
},

{
    "location": "index.html#Table-of-Contents-1",
    "page": "Home",
    "title": "Table of Contents",
    "category": "section",
    "text": "  pages = [\n    \"Home\"       => \"index.md\",\n    \"Overview\"   => \"man/overview.md\",\n    \"Modeling\"   => \"man/modeling.md\",\n    \"Algorithms\" => \"man/algorithms.md\",\n    \"Examples\"   => \"man/examples.md\",\n    \"Benchmarks\" => \"man/benchmarks.md\",\n    \"Developers\" => \"man/developers.md\",\n    \"References\" => \"man/references.md\"\n  ]"
},

{
    "location": "man/overview.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "man/overview.html#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": ""
},

{
    "location": "man/overview.html#Assumptions-1",
    "page": "Overview",
    "title": "Assumptions",
    "category": "section",
    "text": "Given a system of interacting particles and reactions, we assume:the entire particle population is well-mixed (a reaction will happen),\neach reaction depends only on the current system state,\nreaction propensities follow the law of mass action."
},

{
    "location": "man/overview.html#Continuous-Time-Markov-Chains-1",
    "page": "Overview",
    "title": "Continuous-Time Markov Chains",
    "category": "section",
    "text": "The underlying mathematical structure in the present modeling approach is a continuous-time Markov chain mathbfX_t. The random vector mathbfX_t has non-negative components X_ti that count the number of particles of type i at time t. A reaction channel j is characterized by its propensity function r_j(mathbfX_t) that depends on the current vector counts. In a small interval of length Delta t, we expect r_j(mathbfx)Delta t + o(Delta t) reactions of type j to occur. Each reaction j changes the system state by a fixed integer vector mathbfv^j.From the Markov Chain perspective, the chain waits an exponential length of time proportional to the cumulative intensity r_0(mathbfx) = sum_j=1^c r_j(mathbfx). Once a jump occurs, a reaction channel j fires with probability r_j(mathbfx)r_0(mathbfx).Let p_xy(t) denote the finite-time transition probability of going from state x at time 0 to state y at time t. The probability of transition from x to y during a small interval Delta t can be approximated to leading order:\\[ p_{xy}(t + \\Delta t) = p_{xy}(t)\\left[ 1 - \\sum_{j=1}^{c} r_{j}(y) \\Delta t \\right] + \\sum_{j=1}^{c} p_{x,y-v^{j}}(t) r_{j}(y-v^{j}) \\Delta t + o(\\Delta t). \\]One can then form a difference quotient and take the limit Delta t to 0 to arrive at the master equation:\\[ \\frac{\\mathrm{d}}{\\mathrm{d}t} p_{xy}(t) = \\sum_{j=1}^{c} \\left[ p_{x,y-v^{j}}(t) r_{j}(y - v^{j}) - p_{xy}(t) r_{j}(y) \\right], \\]which forms a potentially infinite system of coupled ordinary differential equations. The core of BioSimulator and stochastic simulation in general is to sample the master equation to arrive at statistically exact results."
},

{
    "location": "man/overview.html#Example-Reaction-Channels-1",
    "page": "Overview",
    "title": "Example Reaction Channels",
    "category": "section",
    "text": "Name Reaction r(boldsymbolx) v\nImmigration emptyset to S_1 a_1  1  0  0  0\nDecay S_1 to emptyset a_2x_1 -1 0  0  0\nDimerization S_1 + S_1 to S_2 a_3 binomx_12 -2  1  0  0\nIsomerization S_1 to S_2 a_4 x_1 -1  1  0  0\nDissociation S_2 to S_1 + S_1 a_5 x_2  2 -1  0  0\nBudding S_1 to S_1 + S_1 a_6 x_1  1  0  0  0\nReplacement S_1 + S_2 to S_2 + S_2 a_7 x_1 x_2 -1  1  0  0\nComplex Reaction S_1 + S_2 to S_3 + S_4 a_8 x_1 x_2 -1  -1  1  1"
},

{
    "location": "man/modeling.html#",
    "page": "Modeling",
    "title": "Modeling",
    "category": "page",
    "text": ""
},

{
    "location": "man/modeling.html#Modeling-1",
    "page": "Modeling",
    "title": "Modeling",
    "category": "section",
    "text": "CurrentModule = BioSimulator"
},

{
    "location": "man/modeling.html#Creating-a-Model-1",
    "page": "Modeling",
    "title": "Creating a Model",
    "category": "section",
    "text": "The Network type is the starting point of the modeling process in BioSimulator. It represents a collection related Species that interact through the rules defined by Reactions."
},

{
    "location": "man/modeling.html#BioSimulator.Network",
    "page": "Modeling",
    "title": "BioSimulator.Network",
    "category": "Type",
    "text": "Network(id)\n\nConstruct an empty Network representing a system of interacting particles. A particle is represented by a Species, and an interaction is represented by a Reaction.\n\nAdd a Species or Reaction using <=:\n\nm <= Species(\"X\", 100) # Adds a Species\nm <= Reaction(\"birth\", 2.0, \"X --> X + X\")\n\nSimulate a Network by calling simulate:\n\nsimulate(m, SSA, time=1000.0, epochs=1000, trials=100)\n\nVisualize a Network as a Petri net:\n\nvisualize(m)\n\nSee also: Species, Reaction, simulate, visualize\n\nArguments\n\nid: An string identifier for the Network.\n\n\n\n"
},

{
    "location": "man/modeling.html#BioSimulator.Species",
    "page": "Modeling",
    "title": "BioSimulator.Species",
    "category": "Type",
    "text": "Species(id, [value=0])\n\nDefine a Species with a name and initial population value.\n\n\n\n"
},

{
    "location": "man/modeling.html#BioSimulator.Reaction",
    "page": "Modeling",
    "title": "BioSimulator.Reaction",
    "category": "Type",
    "text": "Reaction(name, rate, formula)\n\nDefine a Reaction with a name, a stochastic rate constant, and a formula.\n\nThe formula should use Species from a Network. Multiple reactants are separated by a + symbol; similarly for products. Reactants and products are separated by a --> symbol. Refer to the examples below:\n\nExamples\n\nReaction(\"birth\", 1.0, \"X –> X + X\")\nReaction(\"death\", 0.5, \"X –> 0\") # no products\nReaction(\"immigration\", 0.5, \"0 –> X\") # no reactants\n\n\n\n"
},

{
    "location": "man/modeling.html#Interface-1",
    "page": "Modeling",
    "title": "Interface",
    "category": "section",
    "text": "  Network  Species  Reaction"
},

{
    "location": "man/modeling.html#BioSimulator.simulate",
    "page": "Modeling",
    "title": "BioSimulator.simulate",
    "category": "Function",
    "text": "simulate{T}(model, [algorithm::Type{T}=SSA];\n  time   :: AbstractFloat=1.0,\n  epochs :: Integer=1,\n  trials :: Integer=1,\n  algvars...)\n\nSimulate a Network using the given algorithm.\n\nThe simulation routine will run until the termination time and record the system state at evenly spaced epochs. This repeats for the given number of trials.\n\nArguments\n\nmodel: The Network to simulate.\nalgorithm: The algorithm used to carry out simulation. One of SSA, FRM, NRM, ODM, or SAL.\n\nOptional Arguments\n\ntime: The amount of time to simulate the model, in units of the model.\nepochs: The number of times to sample the vector of counts.\ntrials: The number of independent realizations to generate.\nkwargs: Additional keyword arguments specific to each algorithm.\n\n\n\n"
},

{
    "location": "man/modeling.html#Running-Simulations-1",
    "page": "Modeling",
    "title": "Running Simulations",
    "category": "section",
    "text": "  simulate"
},

{
    "location": "man/modeling.html#Visualizing-Results-1",
    "page": "Modeling",
    "title": "Visualizing Results",
    "category": "section",
    "text": "Simulation output is stored as a 3D array within a PartialHistory type wherethe i dimension spans the species in the Network,\nthe j dimension spans the epochs sampled, and\nthe k dimension spans the number of trials.BioSimulator provides quick visualizations of sample and mean trajectories, distributions, and phase plots via the Plots package. Please see the Examples page for specific examples."
},

{
    "location": "man/modeling.html#Manipulating-and-Exporting-Data-1",
    "page": "Modeling",
    "title": "Manipulating and Exporting Data",
    "category": "section",
    "text": "Manipulating the simulation data directly can be cumbersome. BioSimulator uses the DataFrames package to extract results as an easy-to-use DataFrame. We note that Plots supports plotting DataFrame objects. One can write simulation results to disk directly, but we recommend using the interface provided by DataFrames."
},

{
    "location": "man/modeling.html#GUI-Generation-1",
    "page": "Modeling",
    "title": "GUI Generation",
    "category": "section",
    "text": "BioSimulator provides experimental automatic graphical user interface generation. One must specify the parameters to expose in the interface and provide default values where necessary."
},

{
    "location": "man/algorithms.html#",
    "page": "Algorithms",
    "title": "Algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "man/algorithms.html#BioSimulator.SSA",
    "page": "Algorithms",
    "title": "BioSimulator.SSA",
    "category": "Type",
    "text": "SSA\n\nGillespie's Direct Method (SSA). Simulates a system of coupled reactions and species by computing the time to the next reaction and searching on the CMF.\n\nInternals\n\nend_time: The termination time, supplied by a user.\nt: The current simulation time.\n\n\n\n"
},

{
    "location": "man/algorithms.html#BioSimulator.FRM",
    "page": "Algorithms",
    "title": "BioSimulator.FRM",
    "category": "Type",
    "text": "FRM\n\nGillespie's First Reaction Method. It is statistically equivalent to SSA, but more computationally expensive: it computes the time to the next reaction as the minimum waiting time relative to the next firing times of each reaction.\n\nInternals\n\nend_time: The termination time, supplied by a user.\nt: The current simulation time.\n\n\n\n"
},

{
    "location": "man/algorithms.html#BioSimulator.NRM",
    "page": "Algorithms",
    "title": "BioSimulator.NRM",
    "category": "Type",
    "text": "NRM\n\nGibson and Bruck's Next Reaction Method, statistically equivalent to SSA. It provides better computational efficiency on networks with loosely connected reactions.\n\nInternals\n\nend_time: The termination time, supplied by a user.\nt: The current simulation time.\npq: A priority queue that sorts reaction according to the their next firing times.\n\n\n\n"
},

{
    "location": "man/algorithms.html#BioSimulator.ODM",
    "page": "Algorithms",
    "title": "BioSimulator.ODM",
    "category": "Type",
    "text": "ODM\n\nOptimized Direct Method. Similar to SSA, with the added benefit of sorting reactions according to their propensities over time. This improves the search on the CMF when selecting the next reaction to fire.\n\nInternals\n\nend_time: The termination time, supplied by a user.\nt: The current simulation time.\n\n\n\n"
},

{
    "location": "man/algorithms.html#BioSimulator.SAL",
    "page": "Algorithms",
    "title": "BioSimulator.SAL",
    "category": "Type",
    "text": "SAL\n\nStep Anticipation-τ Leaping. An algorithm method that improves upon the accuracy of τ-leaping techniques by incorporating a first-order Taylor expansion.\n\nInternals\n\nend_time: The termination time, supplied by a user.\nt: The current simulation time.\nϵ: A parameter controlling the size of leaps. Higher values allow for larger leaps, but may compromise the accuracy of results.\nδ: A τ-leaping parameter used to switch between SSA and SAL. For example, if intensity < δ the algorithm carries out an ordinary SSA step to avoid negative species populations.\nα: A parameter used to contract a τ-leap in the event of negative populations. Aggresive contraction (lower values) will bias sample paths.\ndxdt: Rates of change for each species.\ndrdt: Rates of change for each reaction propensity.\nevents: The number of Poisson arrivals within a τ-leap interval.\n\n\n\n"
},

{
    "location": "man/algorithms.html#Algorithms-1",
    "page": "Algorithms",
    "title": "Algorithms",
    "category": "section",
    "text": "CurrentModule = BioSimulatorSSA\nFRM\nNRM\nODM\nSAL"
},

{
    "location": "man/examples.html#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "man/examples.html#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": ""
},

{
    "location": "man/examples.html#Birth-Death-Immigration-Process-1",
    "page": "Examples",
    "title": "Birth-Death-Immigration Process",
    "category": "section",
    "text": "using BioSimulatorKendall's process is a birth-death-immigration process describing the dynamics of a population using a continuous-time Markov chain. Individuals in the population behave as particles that reproduce at a rate alpha, decay at a rate mu, and immigrate into the population at a rate nu.  model = Network(\"Kendall's Process\")\n\n  model <= Species(\"X\", 5)\n\n  model <= Reaction(\"birth\", 2.0, \"X --> X + X\")\n  model <= Reaction(\"death\", 1.0, \"X --> 0\")\n  model <= Reaction(\"immigration\", 0.5, \"0 --> X\")\n\n  fig = visualize(model)\n  using TikzPictures; save(SVG(\"kendallfig\"), fig) # hide(Image: )  result = simulate(model, algorithm=SSA, time=4.0, epochs=40, trials=1000)\n\n  plot(\n    meantrajectory(result),\n    freqhistogram(result, 4.0),\n    layout = 2,\n    size   = (800,400)\n  )\n\n  savefig(\"kendall.svg\"); nothing # hide(Image: )"
},

{
    "location": "man/examples.html#Enzyme-Kinetics-1",
    "page": "Examples",
    "title": "Enzyme Kinetics",
    "category": "section",
    "text": "using BioSimulatorMichaelis-Menten enzyme kinetics is a stepwise process combining first- and second order reactions to describe the conversion of a substrate into a product. An enzyme E binds to a substrate S to form a complex SE. Conversion does not happen immediately, so SE may revert to its two components or result in a product P and enzyme E.  model = Network(\"enzyme kinetics\")\n\n  model <= Species(\"S\", 301)\n  model <= Species(\"E\", 100)\n  model <= Species(\"SE\",  0)\n  model <= Species(\"P\",   0)\n\n  model <= Reaction(\"Binding\", 0.00166, \"S + E --> SE\")\n  model <= Reaction(\"Dissociation\", 0.0001, \"SE --> S + E\")\n  model <= Reaction(\"Conversion\", 0.1, \"SE --> P + E\")\n\n  fig = visualize(model)\n  using TikzPictures; save(SVG(\"mmekfig\"), fig) # hide(Image: )  result = simulate(model, algorithm=SSA, time=50.0, epochs=100, trials=1000)\n\n  plot(\n    meantrajectory(result),\n    freqhistogram(result, 50.0, select=[\"S\"]),\n    freqhistogram(result, 50.0, select=[\"E\"]),\n    freqhistogram(result, 50.0, select=[\"SE\"]),\n    freqhistogram(result, 50.0, select=[\"P\"]),\n    size   = (800, 600),\n    layout = @layout [a{0.5h}\n                      grid(2,2)]\n  )\n\n  savefig(\"mmek.svg\"); nothing # hide(Image: )"
},

{
    "location": "man/examples.html#Auto-Regulatory-Gene-Network-1",
    "page": "Examples",
    "title": "Auto-Regulatory Gene Network",
    "category": "section",
    "text": "using BioSimulatorThe influence of noise at the cellular level is difficult to capture in deterministic models. Stochastic simulation is appropriate for the study of regulatory mechanisms in genetics, where key species may be present in low numbers.  function autoreg(;k1=1.0, k1r=10.0, k2=0.01, k3=10.0, k4=1.0, k4r=1.0, k5=0.1, k6=0.01)\n    model = Network(\"auto-regulation\")\n\n    model <= Species(\"gene\",   10)\n    model <= Species(\"P2_gene\", 0)\n    model <= Species(\"RNA\",     0)\n    model <= Species(\"P\",       0)\n    model <= Species(\"P2\",      0)\n\n    model <= Reaction(\"repression binding\", k1, \"gene + P2 --> P2_gene\")\n    model <= Reaction(\"reverse repression binding\", k1r, \"P2_gene --> gene + P2\")\n    model <= Reaction(\"transcription\", k2, \"gene --> gene + RNA\")\n    model <= Reaction(\"translation\", k3, \"RNA --> RNA + P\")\n    model <= Reaction(\"dimerization\", k4, \"P + P --> P2\")\n    model <= Reaction(\"dissociation\", k4r, \"P2 --> P + P\")\n    model <= Reaction(\"RNA degradation\", k5, \"RNA --> 0\")\n    model <= Reaction(\"protein degradation\", k6, \"P --> 0\")\n\n    return model\n  end\n\n  model = autoreg()\n\n  fig = visualize(model)\n  using TikzPictures; save(SVG(\"genefig\"), fig) # hide(Image: )  result = simulate(model, algorithm=SSA, time=1000.0, epochs=500, trials=100)\n\n  plot(\n    meantrajectory(result, select=[\"P\", \"P2\"]),\n    freqhistogram(result, 1000.0, select=[\"P\"]),\n    freqhistogram(result, 1000.0, select=[\"P2\"]),\n    layout = @layout [a{0.5h}\n                      grid(1,2)]\n  )\n\n  savefig(\"gene.svg\"); nothing # hide(Image: )"
},

{
    "location": "man/benchmarks.html#",
    "page": "Benchmarks",
    "title": "Benchmarks",
    "category": "page",
    "text": ""
},

{
    "location": "man/benchmarks.html#Benchmarks-1",
    "page": "Benchmarks",
    "title": "Benchmarks",
    "category": "section",
    "text": ""
},

{
    "location": "man/developers.html#",
    "page": "Developers",
    "title": "Developers",
    "category": "page",
    "text": ""
},

{
    "location": "man/developers.html#Developers-1",
    "page": "Developers",
    "title": "Developers",
    "category": "section",
    "text": "CurrentModule = BioSimulator"
},

{
    "location": "man/developers.html#The-simulate-Routine-1",
    "page": "Developers",
    "title": "The simulate Routine",
    "category": "section",
    "text": "The flexibility of the simulate routine allows one to easily implement additional algorithms:function serial_simulate(...)\n\n  init!(algorithm, Xt, r)\n\n  for i in 1:max_trials\n    copy!(Xt, X0)\n    update_all_propensities!(a, r, Xt)\n    reset!(algorithm, a)\n    interval = 1\n\n    while !done(algorithm)\n      interval = update!(output, Xt, get_time(algorithm), interval, i)\n      step!(algorithm, Xt, r)\n    end\n\n    interval = update!(output, Xt, get_time(algorithm), interval, i)\n  end\n\n  return output\nendHere is a breakdown of what happens:init!(...): This performs any initial computations and data structure setup to avoid allocating memory within a loop.\nfor i in 1:max_trials: This loop iterates over independent simulations.\ncopy!(...): Reset the populations to their initial values.\nreset!(...): Reset algorithm variables and data structures.\nwhile !done(...): The main loop in the algorithm.\nupdate!(...): Sample the state vector and record it as needed.\nstep!(...): Execute one step in the algorithm.Thus, one simply needs to implement the init!, reset!, and step! functions for an algorithm."
},

{
    "location": "man/developers.html#BioSimulator.Algorithm",
    "page": "Developers",
    "title": "BioSimulator.Algorithm",
    "category": "Type",
    "text": "An abstract type that species a common interface for all algorithms implemented in BioSimulator.\n\nAll concrete Algorithm types should be immutable and avoid subtyping Algorithm directly; see ExactMethod and TauLeapMethod. The following are normally implemented for Algorithm:\n\nget_time : Return the current simulation time.\nend_time : Return the termination time specified by the user.\ndone : Has the simulation terminated? Defaults to time >= end_time.\ninit! : Initialize an Algorithm prior to simulation. Memory allocation should happen here because it occurs outside loops in the main simulation routine.\nstep! : Carry out a simulation step. This should (1) update simulation time, (2) update the species counts, (3) update propensities, and (4) update any relevant data structures.\n\n\n\n"
},

{
    "location": "man/developers.html#The-Algorithm-Interface-1",
    "page": "Developers",
    "title": "The Algorithm Interface",
    "category": "section",
    "text": "  Algorithm"
},

{
    "location": "man/developers.html#BioSimulator.ExactMethod",
    "page": "Developers",
    "title": "BioSimulator.ExactMethod",
    "category": "Type",
    "text": "An ExactMethod is an abstract type for algorithms that are equivalent to Gillespie's method for stochastic simulation.\n\nIt subtypes the Algorithm type and provides:\n\nreset! : This method defaults to resetting the interal simulation time for an ExactMethod after each trial.\n\n\n\n"
},

{
    "location": "man/developers.html#The-ExactMethod-Interface-1",
    "page": "Developers",
    "title": "The ExactMethod Interface",
    "category": "section",
    "text": "  ExactMethod"
},

{
    "location": "man/developers.html#BioSimulator.TauLeapMethod",
    "page": "Developers",
    "title": "BioSimulator.TauLeapMethod",
    "category": "Type",
    "text": "A TauLeapMethod is an abstract type for algorithms that utilize the τ-leaping strategy to accelerate simulation.\n\nIt subtypes the Algorithm type and provides:\n\nevents : Returns the number of reaction events.\nreset! : This method defaults to resetting the internal simulation time for a TauLeapMethod after each trial.\n\n\n\n"
},

{
    "location": "man/developers.html#The-TauLeapMethod-Interface-1",
    "page": "Developers",
    "title": "The TauLeapMethod Interface",
    "category": "section",
    "text": "  TauLeapMethod"
},

{
    "location": "man/references.html#",
    "page": "References",
    "title": "References",
    "category": "page",
    "text": ""
},

{
    "location": "man/references.html#References-1",
    "page": "References",
    "title": "References",
    "category": "section",
    "text": ""
},

{
    "location": "man/references.html#Papers-1",
    "page": "References",
    "title": "Papers",
    "category": "section",
    "text": "[1] Jeff Bezanson, Stefan Karpinski, Viral B. Shal, and Alan Edelman. Julia: A fast dynamic language for technical computing. arxiv:1209.5145 [cs.PL], September 2012\n[2] Yang Cao and Linda R. Petzold. Accuracy limitations and the measurement of errors in the stochastic simulation of chemically reacting systems. Journal of Computational Physics, 212: 6-24, 2006.\n[3] Yang Cao, Hong Li, and Linda Petzold. Efficient formulation of the stochastic simulation algorithm for chemically reacting systems. Journal of Chemical Physics, 121(9) 2004.\n[4] Yang Cao, Daniel T. Gillespie, and Linda R. Petzold. Efficient step size selection for the tau leaping simulation method. Journal of Chemical Physics, 124(4), 2006.\n[5] Michael A. Gibson and Jehoshua Bruck. E cient exact stochastic simulation of chemical systems with many species and many channels. Journal of Physical Chemistry, 104(9), 1999.\n[6] Daniel T. Gillespie. A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. Journal of Computational Physics, 22(4):403–434, 1976.\n[7] Daniel T. Gillespie. Approximate accelerated stochastic simulation of chemically reacting systems. Journal of Chemical Physics, 115(4), 2001.\n[8] Daniel T. Gillespie and Linda R. Petzold. Improved leap-size selection for accelerated stochastic simulation. Journal of Chemical Physics, 19:8229–8234, 2003.\n[9] Desmond J. Higham. Modeling and simulating chemical reactions. SIAM review, 50(2), 2008.\n[10] Sean Mauch. E cient formulations for exact stochastic simulation of chemical systems. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 8(1):27–35, 2011.\n[11] Hana El Samad, M Khammash, Linda R. Petzold, and Daniel T. Gillespie. Stochastic modeling of gene regulatory networks. International Journal of Robust and Nonlinear Control, 00:1–6, 2002.\n[12] Kevin R. Sanft, Sheng Wu, Min Roh, Jin Fu, Rone Kwei Lim, and Linda R. Petzold. Stochkit2: Software for discrete stochastic simulation of biochemical systems with events. Bioinformatics, 27(17):2457–2458, 2011.\n[13] Mary E. Sehl, Alexander L. Alekseyenko, and Kenneth L. Lange. Accurate stochastic simulation via the step anticipation tau-leaping (SAL) algorithm. Journal of Computational Biology, 16:1195–1208, 2009.\n[14] Timothy G. Vaughan and Alexei J. Drummond. A stochastic simulator of birth-death master equations with applications to phylodynamics. Molecular Biology and Evolution, 30(6):1480– 1493, 2013.\n[15] Darren J. Wilkinson. Stochastic Modeling for Systems Biology. Chapman & Hall/CRC Press, Boca Raton, FL, 1 edition, 2006."
},

{
    "location": "man/references.html#Julia-Packages-1",
    "page": "References",
    "title": "Julia Packages",
    "category": "section",
    "text": "DataFrames\nPlots\nInteract\nIterators\nIJulia\nDistributions\nReexport\nTikzGraphs\nLightGraphs"
},

]}
