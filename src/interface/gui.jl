# species widgets
function make_species_widgets(model, list_of_species)
  species = species_list(model)
  species_controls = Interact.Textbox{Int64}[]

  for key in list_of_species
    s   = species[symbol(key)]
    val = s.population
    box = textbox(value=val, label=key)
    push!(species_controls, box)
  end

  return species_controls
end

function add_species_callback(gui, species_controls)
  for ctrl in species_controls
    id = ctrl.label
    gui = add_species_callback(gui, ctrl, id)
  end
  return gui
end

function add_species_callback(gui, ctrl, id)
  gui = map((m, val) -> m <= Species(id, val), signal(gui), signal(ctrl))
end

# reaction widgets
function make_reaction_widgets(model, list_of_reactions)
  reactions = reaction_list(model)
  reaction_controls = Interact.Textbox{Float64}[]

  for key in list_of_reactions
    r   = reactions[symbol(key)]
    val = r.rate
    box = textbox(value=val, label="$key")
    push!(reaction_controls, box)
  end

  return reaction_controls
end

function add_reaction_callback(gui, reaction_controls)
  for ctrl in reaction_controls
    r    = reaction_list(value(gui))[symbol(ctrl.label)]
    id   = r.id
    reac = r.reactants
    prod = r.products

    gui = add_reaction_callback(gui, ctrl, id, reac, prod)
  end
  return gui
end

function add_reaction_callback(gui, ctrl, id, r, p)
  gui = map((m, val) -> m <= Reaction(id, val, r, p), signal(gui), signal(ctrl))
end

function init_and_simulate(vars)
  model  = vars[1]
  algnm  = vars[2]
  time   = vars[3]
  trials = vars[4]
  epochs = vars[5]

  algorithm = algnm(end_time=time)

  return simulate(model, algorithm, sampling_interval=time/epochs, nrlz=trials)
end

function make_visualization(fplot::Function, result, selected)
  return fplot(result, select=[selected])
end

function generate_gui(model::Network, list_of_species, list_of_reactions)
  # make the Network into a signal
  m_sig = Signal(model)

  # make species controls
  swidgs = make_species_widgets(value(m_sig), list_of_species)

  # make reaction rate controls
  rwidgs = make_reaction_widgets(value(m_sig), list_of_reactions)

  # make algorithm controls
  runbtn = button("Run")

  sel_alg    = dropdown(ALGORITHMS)
  set_time   = textbox(value=1.0, label="time")
  set_trials = textbox(value=1,   label="trials")
  set_epochs = textbox(value=1,   label="epochs")

  # display everything so far
  map(display, swidgs)
  map(display, rwidgs)
  display(sel_alg)
  display(set_time)
  display(set_trials)
  display(set_epochs)
  display(runbtn)

  add_species_callback(m_sig, swidgs)
  add_reaction_callback(m_sig, rwidgs)

  # couple updates with a button click
  vars = map((m, alg, t, k, n) -> (m, alg, t, k, n),
  signal(m_sig), signal(sel_alg), signal(set_time), signal(set_trials), signal(set_epochs))

  update_sig = sampleon(signal(runbtn), vars)
  result = map(init_and_simulate, signal(update_sig))

  # make plot type control
  sel_plot = dropdown(PLOT_TYPES)

  # make visualization control
  sel_species = dropdown(collect(keys(value(result).id2ind)))
  sel_epoch      = textbox(value=1, label="epoch no.")

  display(sel_plot)
  display(sel_species)

  # couple result to visualization
  viz = map(make_visualization, signal(sel_plot), signal(result), signal(sel_species))

  return viz
end
