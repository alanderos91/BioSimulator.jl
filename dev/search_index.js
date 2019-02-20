var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#BioSimulator.jl-1",
    "page": "Home",
    "title": "BioSimulator.jl",
    "category": "section",
    "text": "A stochastic simulation framework for Julia."
},

{
    "location": "#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "Many complex systems in biology are analytically intractable, and dynamical predictions based on deterministic models can be grossly misleading. Stochastic simulation algorithms based on continuous-time Markov chains allow researchers to generate accurate time-evolution trajectories, test the sensitivity of models to key parameters, and quantify frequencies of rare events.Situations where stochastic simulation is especially helpful involve:rare events such as extinction and mutation,\nkey molecules present in small numbers,\nrare reactions with dramatic influence on the dynamics of the system\npopulation cycles arising from demographic stochasticity.Examples of such systems include gene expression networks, tumor suppressor pathways, and demographic and ecological systems.BioSimulator.jl aims to provide researchers interested in such phenomena with a fast, reliable, user-friendly, and open-source modeling tool in Julia."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "BioSimulator.jl must be installed with Pkg.clone in the Julia REPL:Pkg.clone(\"https://github.com/alanderos91/biosimulator.jl.git\", \"BioSimulator\")You can start using BioSimulator.jl in scripts or the REPL with the command:using BioSimulator"
},

{
    "location": "#Additional-tools-1",
    "page": "Home",
    "title": "Additional tools",
    "category": "section",
    "text": ""
},

{
    "location": "#DataFrames-1",
    "page": "Home",
    "title": "DataFrames",
    "category": "section",
    "text": "The SimulationSummary returned by simulate can be converted into a DataFrame using the DataFrames.jl package. The conversion is straightforward after installing the package with the command Pkg.add(\"DataFrames\"):# make sure the package is loaded\nusing DataFrames\n\n# simulate a model and save the results\nresult = simulate(model...)\n\n# returns a DataFrame\nDataFrame(result)The first two columns indicate the time point and trial of a record (row). The remaining columns represent the species counts and are labeled using the given names."
},

{
    "location": "#Plotting-1",
    "page": "Home",
    "title": "Plotting",
    "category": "section",
    "text": "The plotting defaults provided by BioSimulator.jl require the Plots.jl package. You can install it withPkg.add(\"Plots\")BioSimulator.jl does not load the Plots.jl package by default. Any time you need plotting functionality, simply load the package:# if BioSimulator is already loaded\nusing Plots\n\n# if you\'re just starting\nusing BioSimulator, PlotsNote that Plots.jl is independent of BioSimulator.jl and can be used without BioSimulator.jl. Please consult the Plots.jl documentation for additional details."
},

{
    "location": "#Petri-nets-1",
    "page": "Home",
    "title": "Petri nets",
    "category": "section",
    "text": "BioSimulator.jl should install the TikzGraphs.jl package by default. You can try generating a Petri net in a Jupyter notebook for a model using visualize(model). If you want to generate the figure in a script and save it to a file, install the TikzPictures.jl package. You can then save a figure using:import TikzPictures: save\n\nfigure = visualize(model)\nsave(PDF(filename), figure)"
},

{
    "location": "#Table-of-Contents-1",
    "page": "Home",
    "title": "Table of Contents",
    "category": "section",
    "text": "  pages = [\n    \"Home\"       => \"index.md\",\n    \"Overview\"   => \"man/overview.md\",\n    \"Algorithms\" => \"man/algorithms.md\",\n    \"Examples\"   => \"man/examples.md\",\n  ]"
},

{
    "location": "man/overview/#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "man/overview/#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": "CurrentModule = BioSimulator"
},

{
    "location": "man/overview/#Creating-a-Model-1",
    "page": "Overview",
    "title": "Creating a Model",
    "category": "section",
    "text": "The Network type is the starting point of the modeling process in BioSimulator.jl. A Network is a system of coupled Species and Reaction channels."
},

{
    "location": "man/overview/#BioSimulator.Network",
    "page": "Overview",
    "title": "BioSimulator.Network",
    "category": "type",
    "text": "Network(id)\n\nConstruct an empty Network representing a system of interacting particles. A particle is represented by a Species, and an interaction is represented by a Reaction.\n\nAdd a Species or Reaction using <=:\n\nm <= Species(\"X\", 100) # Adds a Species\nm <= Reaction(\"birth\", 2.0, \"X --> X + X\")\n\n\n\n\n\n"
},

{
    "location": "man/overview/#BioSimulator.Species",
    "page": "Overview",
    "title": "BioSimulator.Species",
    "category": "type",
    "text": "Species(id, [value=0])\n\nDefine a Species with a name and initial population value.\n\n\n\n\n\n"
},

{
    "location": "man/overview/#BioSimulator.Reaction",
    "page": "Overview",
    "title": "BioSimulator.Reaction",
    "category": "type",
    "text": "Reaction(name, rate, formula)\n\nDefine a Reaction with a name, a stochastic rate constant, and a formula.\n\nThe formula should use Species from a Network. Multiple reactants are separated by a + symbol; similarly for products. Reactants and products are separated by a --> symbol. Refer to the examples below:\n\nExamples\n\nReaction(\"birth\", 1.0, \"X –> X + X\")\nReaction(\"death\", 0.5, \"X –> 0\") # no products\nReaction(\"immigration\", 0.5, \"0 –> X\") # no reactants\n\n\n\n\n\n"
},

{
    "location": "man/overview/#Interface-1",
    "page": "Overview",
    "title": "Interface",
    "category": "section",
    "text": "  Network  Species  Reaction"
},

{
    "location": "man/overview/#BioSimulator.simulate",
    "page": "Overview",
    "title": "BioSimulator.simulate",
    "category": "function",
    "text": "simulate(model::Network, algname::T, [output_type::Val{S} = Val(:fixed)]; time::Float64=1.0, epochs::Int=1, trials::Int=1, track_stats::Bool=false, kwargs...) where {T,S}\n\nSimulate a model with algname. The simulation routine will run until the termination time for the given number of trials.\n\nArguments\n\nmodel: The Network to simulate.\nalgname: The name of the algorithm to carry out the simulation. One of Direct(), FirstReaction(), NextReaction(), OptimizedDirect(), TauLeaping(), or StepAnticipation().\noutput_type = Val(:fixed): The type of output to record. The Val(:full) option records every simulation step, whereas the Val(:fixed) option records a fixed number of epochs.\n\nOptional Arguments\n\ntime = 1.0: The amount of time to simulate the model, in units of the model.\nepochs = 1: The number of times to sample the vector of counts.\ntrials = 1: The number of independent realizations to generate.\ntrack_stats = false: An option to toggle tracking algorithm statistics (e.g. number of steps).\nkwargs: Additional keyword arguments specific to each algorithm. Reference a specific algname for details.\n\n\n\n\n\n"
},

{
    "location": "man/overview/#Running-Simulations-1",
    "page": "Overview",
    "title": "Running Simulations",
    "category": "section",
    "text": "  simulate"
},

{
    "location": "man/overview/#Parallelization-1",
    "page": "Overview",
    "title": "Parallelization",
    "category": "section",
    "text": "One can start Julia with the --procs=N option to enable parallelization. Here N is the number of threads available. In this case, the simulate routine will automatically delegate individual simulations to each available thread until the total number of trials is satisfied."
},

{
    "location": "man/overview/#Simulation-results-1",
    "page": "Overview",
    "title": "Simulation results",
    "category": "section",
    "text": "The result of simulate is a SimulationSummary that stores information about a model, the algorithm used in a simulation, and simulation results. We summarize the fields of a result of type SimulationSummary below:result.model: The Network that was simulated.\nresult.algorithm_name: The name of the algorithm used in the simulation.\nresult.algorithm_params: The settings for the algorithm, including the epochs, time span, and trials used.\nresult.algorithm_stats: Any statistics for the algorithm used. Direct methods simply record the number of steps taken over all trials. Tau-leaping methods also record the number of negative excursions and the number of times the recovery mechanism is used.\nresult.simulation_data: The results of the simulation. It is a vector whose elements are SamplePaths (using the Val(:full) option) or RegularPaths (using the Val(:fixed) option).\nresult.id2index: A dictionary that maps the name of a species to an index into the state vector mathbfX_t."
},

{
    "location": "man/overview/#SamplePath-output-1",
    "page": "Overview",
    "title": "SamplePath output",
    "category": "section",
    "text": "Using the Val(:full) option for simulate records the state vector of the process after every event. In this case, the simulation_data field is a vector of SamplePaths. Each entry of the vector is an individual realization of the process. A SamplePath is a variable-length collection of records at various time points. Given a SamplePath, say xw, one can access the time points by xw.tdata and the states by xw.xdata. Specifically, the i-th record corresponds to a state mathbfX_t_i at time t_i.The advantage of SamplePath output and the Val(:full) option is that one can recover the distribution of each species at different points in time."
},

{
    "location": "man/overview/#RegularPath-output-1",
    "page": "Overview",
    "title": "RegularPath output",
    "category": "section",
    "text": "Using the Val(:fixed) option for simulate records the state vector at pre-defined intervals determined by the epochs option. Given that epochs = n, the time span of a simulation is broken up into n intervals of equal length. After each step, the simulation checks to see if it has stepped over into the next epoch. If it has, then it records the last state as the value of the stochastic process at the previous epoch. This option produces a vector of RegularPaths, with each entry corresponding to an independent realization of the process. Each RegularPath has the fields tdata and xdata representing time points and states, respectively. Importantly, xdata is a pre-allocated matrix with column i storing the value mathbfX_t_i.One advantage of RegularPath output and the Val(:fixed) option is a slight performance boost attributed to minimal dynamic memory allocation. This advantage may disappear if one uses far more epochs than the expected number of events in a given time span. However, care must be taken in choosing the number of epochs. Using too few epochs may grossly simplify the observed qualitative behavior in visualizations. Moreover, only the last save point may be used to reliably estimate a probability distribution. Another advantage of this option is that RegularPaths facilitate computing summary statistics, including mean trajectories."
},

{
    "location": "man/overview/#Plotting-1",
    "page": "Overview",
    "title": "Plotting",
    "category": "section",
    "text": "BioSimulator.jl provides plot recipes that operates on a SimulationSummary. The basic interface is as follows:plot(result::SimulationSummary{T}; plot_type = :trajectory,\n    trial = nothing, species = nothing, epochs = nothing) where TThere are three different options for plot_type: :trajectory, :meantrajectory, and :histogram. Each of these options may be modulated by the trial, species, and epochs options."
},

{
    "location": "man/overview/#Trajectory-1",
    "page": "Overview",
    "title": "Trajectory",
    "category": "section",
    "text": "The plot_type = :trajectory option produces the time evolution of a single realization of the underlying stochastic process. The trial option can be used to select different realizations. The species option can be used to specify a vector of species to plot in a figure."
},

{
    "location": "man/overview/#Mean-Trajectory-1",
    "page": "Overview",
    "title": "Mean Trajectory",
    "category": "section",
    "text": "The plot_type = :meantrajectory option produces the averaged time evolution of the underlying stochastic process. The expected value of a particular species is plotted as a point at evenly spaced intervals. Error bars denote one standard deviation away from the mean. Note that this option is fully compatible with the Val(:full) option by discretizing the time span. This is achieved by passing a value for the epochs option. The species option can be used to specify a vector of species to plot in a figure."
},

{
    "location": "man/overview/#Histogram-1",
    "page": "Overview",
    "title": "Histogram",
    "category": "section",
    "text": "The plot_type = :histogram option produces an estimate of the distribution of the process at the final time point in the simulation. The species option can be used to specify a vector of species to plot in a figure."
},

{
    "location": "man/algorithms/#",
    "page": "Algorithms",
    "title": "Algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "man/algorithms/#BioSimulator.Direct",
    "page": "Algorithms",
    "title": "BioSimulator.Direct",
    "category": "type",
    "text": "Direct()\n\nGillespie\'s Direct Method. Simulates a system of coupled reactions and species by computing the time to the next reaction and searching on the CMF.\n\nReferences\n\nDaniel T Gillespie. A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. Journal of Computational Physics, 1976. https://doi.org/10.1016/0021-9991(76)90041-3\n\nDaniel T. Gillespie. Exact stochastic simulation of coupled chemical reactions. The Journal of Physical Chemistry, 1977. https://doi.org/10.1021/j100540a008\n\n\n\n\n\n"
},

{
    "location": "man/algorithms/#BioSimulator.FirstReaction",
    "page": "Algorithms",
    "title": "BioSimulator.FirstReaction",
    "category": "type",
    "text": "FirstReaction()\n\nGillespie\'s First Reaction Method.\n\nStatistically equivalent to SSA, but more computationally expensive: it computes the time to the next reaction as the minimum waiting time relative to the next firing times of each reaction.\n\nReferences\n\nDaniel T Gillespie. A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. Journal of Computational Physics, 1976. https://doi.org/10.1016/0021-9991(76)90041-3\n\n\n\n\n\n"
},

{
    "location": "man/algorithms/#BioSimulator.NextReaction",
    "page": "Algorithms",
    "title": "BioSimulator.NextReaction",
    "category": "type",
    "text": "NextReaction()\n\nGibson and Bruck\'s Next Reaction Method. It provides better computational efficiency on networks with loosely connected reactions.\n\nReferences\n\nMichael A. Gibson and Jehoshua Bruck. Efficient exact stochastic simulation of chemical systems with many species and many channels. Journal of Physical Chemistry, 1999. https://doi.org/10.1021/jp993732q\n\n\n\n\n\n"
},

{
    "location": "man/algorithms/#BioSimulator.OptimizedDirect",
    "page": "Algorithms",
    "title": "BioSimulator.OptimizedDirect",
    "category": "type",
    "text": "OptimizedDirect()\n\nOptimized Direct Method. Similar to SSA, with the added benefit of sorting reactions according to their propensities over time. This improves the search on the CMF. Sorting is isdone with pre-simulation.\n\nReferences\n\nYang Cao, Hong Li, and Linda Petzold. Efficient formulation of the stochastic simulation algorithm for chemically reacting systems. Journal of Chemical Physics, 2004. https://doi.org/10.1063/1.1778376\n\n\n\n\n\n"
},

{
    "location": "man/algorithms/#BioSimulator.TauLeaping",
    "page": "Algorithms",
    "title": "BioSimulator.TauLeaping",
    "category": "type",
    "text": "TauLeaping()\n\nPoisson tau-leaping. This version implements the description by Gillespie and Petzold (2003) and provides safeguards against negative populations.\n\nOptional Arguments\n\nepsilon = 0.03: A parameter controlling the size of leaps. Higher values allow for larger leaps, but may compromise the accuracy of results.\ndelta = 2.0: A parameter used to switch between a Gillespie update and tau-leaping. If the next tau leap is less than some multiple delta of the expected time for one event to occur, then StepAnticipation carries out a Gillespie step to avoid negative species populations and save on efficiency.\nbeta = 0.75: A parameter between 0 and 1 used to contract a tau leap in the event of negative populations. Aggresive contraction (lower values) will bias sample paths.\n\nReferences\n\nDaniel T. Gillespie and Linda R. Petzold. Improved leap-size selection for accelerated stochastic simulation. Journal of Chemical Physics, 2003. https://doi.org/10.1063/1.1613254\n\n\n\n\n\n"
},

{
    "location": "man/algorithms/#BioSimulator.StepAnticipation",
    "page": "Algorithms",
    "title": "BioSimulator.StepAnticipation",
    "category": "type",
    "text": "StepAnticipation()\n\nStep Anticipation tau-Leaping. An algorithm that improves upon the accuracy of Poisson tau-leaping techniques by incorporating a first-order Taylor expansion. It provides safeguards against negative populations.\n\nOptional Arguments\n\nepsilon = 0.03: A parameter controlling the size of leaps. Higher values allow for larger leaps, but may compromise the accuracy of results.\ndelta = 2.0: A parameter used to switch between a Gillespie update and tau-leaping. If the next tau leap is less than some multiple delta of the expected time for one event to occur, then StepAnticipation carries out a Gillespie step to avoid negative species populations and save on efficiency.\nbeta = 0.75: A parameter between 0 and 1 used to contract a tau leap in the event of negative populations. Aggresive contraction (lower values) will bias sample paths.\n\nReferences\n\nMary E. Sehl, Alexander L. Alekseyenko, and Kenneth L. Lange. Accurate stochastic simulation via the step anticipation tau-leaping (SAL) algorithm. Journal of Computational Biology, 2009. https://dx.doi.org/10.1089/cmb.2008.0249\n\n\n\n\n\n"
},

{
    "location": "man/algorithms/#Algorithms-1",
    "page": "Algorithms",
    "title": "Algorithms",
    "category": "section",
    "text": "CurrentModule = BioSimulatorDirect\nFirstReaction\nNextReaction\nOptimizedDirect\nTauLeaping\nStepAnticipation"
},

{
    "location": "man/examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "man/examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "The following examples illustrate BioSimulator.jl\'s interface and features. Each code block assumes BioSimulator.jl and Plots.jl are loaded; that is, using BioSimulator, Plots.using BioSimulator, Plots, TikzPictures\ngr(fmt = :png, dpi = 300)using BioSimulator, Plots, TikzPictures\ngr(fmt = :png, dpi = 300)using BioSimulator, Plots, TikzPictures\ngr(fmt = :png, dpi = 300)using BioSimulator, Plots\ngr(fmt = :png, dpi = 300)"
},

{
    "location": "man/examples/#Birth-Death-Immigration-Process-1",
    "page": "Examples",
    "title": "Birth-Death-Immigration Process",
    "category": "section",
    "text": "Kendall\'s process is a birth-death-immigration process describing the dynamics of a population using a continuous-time Markov chain. Individuals in the population behave as particles that reproduce at a rate alpha, decay at a rate mu, and immigrate into the population at a rate nu."
},

{
    "location": "man/examples/#Model-Definition-1",
    "page": "Examples",
    "title": "Model Definition",
    "category": "section",
    "text": "# initialize\nmodel = Network(\"Kendall\'s Process\")\n\n# species definitions\nmodel <= Species(\"X\", 5)\n\n# reaction definitions\nmodel <= Reaction(\"birth\", 2.0, \"X --> X + X\")\nmodel <= Reaction(\"death\", 1.0, \"X --> 0\")\nmodel <= Reaction(\"immigration\", 0.5, \"0 --> X\")\n\n# Petri net; the \"scale = 2\" argument is optional and can be omitted\n# it is used to pass additional options for the underlying Tikz document\nfig = visualize(model, \"scale = 2\")\nTikzPictures.save(SVG(\"kendall_petri.svg\"), fig) # hide(Image: )"
},

{
    "location": "man/examples/#Sample-Output-1",
    "page": "Examples",
    "title": "Sample Output",
    "category": "section",
    "text": "# simulate 1000 realizations with Gillespie algorithm over [0, 4]\n# state is saved over 40 equal-length intervals (fixed-interval output)\nresult = simulate(model, Direct(), time = 4.0, epochs = 40, trials = 1000)\n\n# first panel - mean trajectory\npanel1 = plot(result, plot_type = :meantrajectory)\n\n# second panel - distribution at t = 4.0\npanel2 = plot(result, plot_type = :histogram)\n\n# combine both panels into a single figure; omit the legend\nplot(panel1, panel2, layout = grid(2, 1), legend = nothing)\n\nsavefig(\"kendall.png\"); nothing # hide(Image: )"
},

{
    "location": "man/examples/#Enzyme-Kinetics-1",
    "page": "Examples",
    "title": "Enzyme Kinetics",
    "category": "section",
    "text": "Michaelis-Menten enzyme kinetics is a stepwise process combining first- and second order reactions to describe the conversion of a substrate into a product. An enzyme E binds to a substrate S to form a complex SE. Conversion does not happen immediately, so SE may revert to its two components or result in a product P and enzyme E."
},

{
    "location": "man/examples/#Model-Definition-2",
    "page": "Examples",
    "title": "Model Definition",
    "category": "section",
    "text": "# initialize\nmodel = Network(\"enzyme kinetics\")\n\n# species definitions\nmodel <= Species(\"S\", 301)\nmodel <= Species(\"E\", 100)\nmodel <= Species(\"SE\",  0)\nmodel <= Species(\"P\",   0)\n\n# reaction definitions\nmodel <= Reaction(\"Binding\", 0.00166, \"S + E --> SE\")\nmodel <= Reaction(\"Dissociation\", 0.0001, \"SE --> S + E\")\nmodel <= Reaction(\"Conversion\", 0.1, \"SE --> P + E\")\n\n# Petri net\nfig = visualize(model)\nTikzPictures.save(SVG(\"mmek_petri.svg\"), fig) # hide(Image: )"
},

{
    "location": "man/examples/#Sample-Output-2",
    "page": "Examples",
    "title": "Sample Output",
    "category": "section",
    "text": "Here we plot the mean trajectory of all 4 species over 0 50 and include their distributions at t = 50. Like the previous example, we\'ll put each figure into separate panels. To match colors between two panels, we\'ll extract the colors from the default palette and set the colors explicitly. By default, every species will appear in a plot. We can select which species should appear using the species = list option, where list is an array of species names. In addition, this option can be used to order the species; here we make use of this feature to impose the order S, E, SE, and then P. This makes it so that the first color in our palette corresponds to S, the second to E, and so on.result = simulate(model, Direct(), time = 50.0, epochs = 100, trials = 1000)\n\n# grab colors from default palette\nmycolors = palette(:default)\n\n# set the order and apply colors\nplot(result, plot_type = :meantrajectory, species = [\"S\", \"E\", \"SE\", \"P\"], palette = mycolors[1:4])\n\nsavefig(\"mmek1.png\"); nothing # hide(Image: )We split the second panel into 4 subfigures in order to assign colors explicitly.p1 = plot(result, plot_type = :histogram, species = [\"S\"],  color = mycolors[1])\np2 = plot(result, plot_type = :histogram, species = [\"E\"],  color = mycolors[2])\np3 = plot(result, plot_type = :histogram, species = [\"SE\"], color = mycolors[3])\np4 = plot(result, plot_type = :histogram, species = [\"P\"],  color = mycolors[4])\n\n# combine the subfigures, omit the default title and legend\nplot(p1, p2, p3, p4, layout = grid(2, 2), legend = nothing, title = \"\")\nsavefig(\"mmek2.png\"); nothing # hide(Image: )Alternatively, if color is not an issue, we can use:plot(result, plot_type = :histogram, layout = grid(2, 2), title = \"\")\nsavefig(\"mmek3.png\"); nothing # hide(Image: )"
},

{
    "location": "man/examples/#Auto-Regulatory-Gene-Network-1",
    "page": "Examples",
    "title": "Auto-Regulatory Gene Network",
    "category": "section",
    "text": "The influence of noise at the cellular level is difficult to capture in deterministic models. Stochastic simulation is appropriate for the study of regulatory mechanisms in genetics, where key species may be present in low numbers."
},

{
    "location": "man/examples/#Model-Definition-3",
    "page": "Examples",
    "title": "Model Definition",
    "category": "section",
    "text": "We can wrap the model definition into a function, called autoreg, that can be called to build the model with different sets of parameters.function autoreg(;k1=1.0, k1r=10.0, k2=0.01, k3=10.0, k4=1.0, k4r=1.0, k5=0.1, k6=0.01)\n  # initialize\n  model = Network(\"auto-regulation\")\n\n  # species definitions\n  model <= Species(\"gene\",   10)\n  model <= Species(\"P2_gene\", 0)\n  model <= Species(\"RNA\",     0)\n  model <= Species(\"P\",       0)\n  model <= Species(\"P2\",      0)\n\n  # reaction definitions\n  model <= Reaction(\"repression binding\", k1, \"gene + P2 --> P2_gene\")\n  model <= Reaction(\"reverse repression binding\", k1r, \"P2_gene --> gene + P2\")\n  model <= Reaction(\"transcription\", k2, \"gene --> gene + RNA\")\n  model <= Reaction(\"translation\", k3, \"RNA --> RNA + P\")\n  model <= Reaction(\"dimerization\", k4, \"P + P --> P2\")\n  model <= Reaction(\"dissociation\", k4r, \"P2 --> P + P\")\n  model <= Reaction(\"RNA degradation\", k5, \"RNA --> 0\")\n  model <= Reaction(\"protein degradation\", k6, \"P --> 0\")\n\n  return model\nend\n\n# build with the default parameter values\nmodel = autoreg()\n\n# Petri net\nfig = visualize(model)\nTikzPictures.save(SVG(\"gene_petri.svg\"), fig) # hide(Image: )"
},

{
    "location": "man/examples/#Sample-Output-3",
    "page": "Examples",
    "title": "Sample Output",
    "category": "section",
    "text": "By default, BioSimulator.jl uses fixed-interval output. We can request that the state vector be saved after every step using the Val(:full) option immediately after specifying the algorithm. This option is compatible with mean trajectories; simply specify epochs = n where n is the number of epochs to use in estimating the required means.# simulate with the full output option\nresult = simulate(model, Direct(), Val(:full), time = 100.0, trials = 100)\n\n# plot the mean trajectory using 25 epochs\nplot(result, plot_type = :meantrajectory, species = [\"P\", \"P2\"], epochs = 25)\n\nsavefig(\"gene1.png\"); nothing # hide(Image: )The distributions at t = 4 for P and P2:plot(result, plot_type = :histogram, species = [\"P\", \"P2\"], layout = (2, 1))\n\nsavefig(\"gene2.png\"); nothing # hide(Image: )"
},

{
    "location": "man/examples/#Brusselator-Cascade-1",
    "page": "Examples",
    "title": "Brusselator Cascade",
    "category": "section",
    "text": "The Brusselator is a theoretical model used to study autocatalytic reactions. The reactions in the model areA to X2X + Y to 3XX + B to Y + DX to Ewhere A, B, and D are chemical species assumed to be constant in concentration; only X and Y vary over time.. The species A and B act as inputs to the synthesis and conversion of X, whereas D and E are byproducts with no impact on the system. The Y species acts as a catalyst to the synthesis of X so that X is autocatalytic. Note that the last reaction can be thought of as the decay of X.One can study the role of stochasticity in chemical reaction cascades by coupling Brusselators. A cascade with N steps is modeled asA to X_12 X_n + Y_n to 3 X_n n = 1ldotsNX_n + B to Y_n + D n = 1ldotsNX_n to X_n+1 n = 1ldotsN-1X_N to Ewhere the X_n to X_n+1 reactions effectively couple each Brusselator. The fixed point of the system is given by the concentrations of A and B: (X_n Y_n) = (A A  B) for each n = 1ldotsN, which is stable for B  A^2 + 1 and unstable for B  A^2 + 1. A limit cycle exists in the unstable case which propagates noise down the steps in the cascade. A deterministic model predicts asymptotic stability at the end of the cascade, but small fluctuations from the fixed point are amplified in the stochastic setting.This example walks through an implementation of the Brusselator cascade, including conversion of deterministic rates to stochastic rates."
},

{
    "location": "man/examples/#Model-Definition-4",
    "page": "Examples",
    "title": "Model Definition",
    "category": "section",
    "text": "model = Network(\"Brusselator\")\n\nN = 20    # number of Brusselators\nV = 100.0 # system volume\n\n# ===== \"Deterministic\" Rates =====\nk1 = 1.0  # buffer rate\nk2 = 1.0  # transition/decay rate\nk3 = 1.0  # conversion rate\nk4 = 1.0  # autocatalytic rate\n\n# ===== \"Stochastic\" Rates =====\n# To model a constant buffer X0 we add a zero-order reaction (like immigration)\n# The stochastic rates have to take into account the system volume\n\nγ1 = k1                    # buffer rate\nγ2 = k2 / V                # transition/decay rate\nγ3 = k3 / V                # conversion rate\nγ4 = 2 * k4 / (V * V * V)  # autocatalytic rate\n\nfor i = 1:N\n  # species definitions\n  model <= Species(\"X$(i)\", 0)\n  model <= Species(\"Y$(i)\", 0)\n\n  # autocatalytic reactions\n  model <= Reaction(\"conversion$(i)\",    γ3, \"X$(i) --> Y$(i)\")\n  model <= Reaction(\"autocatalysis$(i)\", γ4, \"X$(i) + X$(i) + Y$(i) --> X$(i) + X$(i) + X$(i)\")\nend\n\nfor i = 2:N\n  # cascades\n  model <= Reaction(\"cascade$(i)\", γ2, \"X$(i-1) --> X$(i)\")\nend\n\nmodel <= Reaction(\"buffer\", γ1, \"0 --> X1\")\nmodel <= Reaction(\"decay\",  γ2, \"X$(N) --> 0\")\nnothing # hide"
},

{
    "location": "man/examples/#Sample-Output-4",
    "page": "Examples",
    "title": "Sample Output",
    "category": "section",
    "text": "Here we plot the sample paths for X and Y at the first and twentieth stages in the cascade; the species are labeled X_1, Y_1, X_20, and Y_20.result = simulate(model, Direct(), Val(:full), time = 10_000.0)\n\np1 = plot(result, plot_type = :trajectory, species = [\"X1\", \"Y1\"], title = \"brusselator 1\")\n\np2 = plot(result, plot_type = :trajectory, species = [\"X20\", \"Y20\"], title = \"brusselator 20\")\n\nplot(p1, p2, layout = grid(2, 1), legend = :best)\n\nsavefig(\"brusselator.png\"); nothing # hide(Image: )"
},

]}
