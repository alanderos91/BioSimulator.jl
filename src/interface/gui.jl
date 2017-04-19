const ALGORITHMS = [SSA, FRM, NRM, ODM, SAL]
const PLOT_TYPES = ["Mean Trajectory", "Histogram"]

function make_species_widgets(model, list_of_species)
  species = species_list(model)
  species_controls = Interact.Textbox{Int64}[]

  for key in list_of_species
    s   = species[Symbol(key)]
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
    r   = reactions[Symbol(key)]
    val = r.rate
    box = textbox(value=val, label="$key")
    push!(reaction_controls, box)
  end

  return reaction_controls
end

function add_reaction_callback(gui, reaction_controls)
  for ctrl in reaction_controls
    r    = reaction_list(value(gui))[Symbol(ctrl.label)]
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
  algo   = vars[2]
  time   = vars[3]
  trials = vars[4]
  epochs = vars[5]

  return simulate(model, algo, time=time, epochs=epochs, trials=trials)
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

    sel_alg    = dropdown(ALGORITHMS, label="algorithm")
    set_time   = textbox(value=1.0,   label="time")
    set_trials = textbox(value=1,     label="trials")
    set_epochs = textbox(value=1,     label="epochs")

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

    return result
end

function plot_interface(result)
    # make visualization controls
    sel_plot    = dropdown(PLOT_TYPES, label="Plot:")
    sel_species = dropdown(collect(keys(value(result).id2ind)), label="Selected Species:")
    sel_epoch   = dropdown(collect(value(result).t_index), label="Selected Time:")
    plt_btn     = button("Plot")

    display(sel_plot)
    display(sel_species)
    display(sel_epoch)
    display(plt_btn)

    # couple result to visualization
    vars = map( (r, t, s, p) -> (r, t, s, p),
        signal(result),
        signal(sel_epoch),
        signal(sel_species),
        signal(sel_plot)
    )
    update_sig = sampleon(signal(plt_btn), vars)
    viz = map(make_visualization, signal(update_sig))

    return viz
end

function make_visualization(vars)
    result   = vars[1]
    time     = vars[2]
    selected = vars[3]
    ftype    = vars[4]

    if ftype == "Mean Trajectory"
        plot(MeanTrajectory(result, selected))
    elseif ftype == "Histogram"
        plot(Histogram(result, selected, time))
    end
end
