abstract Coupling

immutable Tight <: Coupling end
immutable Loose <: Coupling end

export Tight, Loose

function ossa(model::Simulation, t_final::Float64, output::OutputType, dt::Float64, itr::Int, c::Coupling, steps::Int, samples::Int)
	# Unpack model
	init = model.initial
	spcs = model.state
	rxns = model.rxns
	params = model.param

	n = round(Int, t_final / dt) + 1

	job = SimulationJob(itr)

	g = typeof(c) <: Loose ? init_dep_graph(rxns) : DiGraph()

	# Presimulate to sort reactions according to multiscale property.
	# This will modify spcs and rxns
	init_ossa!(spcs, rxns, params, init, steps, samples)

	for i = 1:itr
		reset!(spcs, init)

		result = init_sr(output, spcs, n)
		ssa_steps = 0

		t = 0.0
		t_next = 0.0
		j = 1

		# Compute propensities
		intensity = compute_propensities!(rxns, spcs, params)

		while t < t_final
			t_next, j = update!(output, result, n, t, t_next, dt, j, spcs)
			τ, intensity = ossa_update!(c, spcs, rxns, params, intensity, g)
			t = t + τ
			ssa_steps = ssa_steps + 1
		end
		update!(output, result, n, t_next, dt, j, spcs)
		job[i] = result
	end
	return job
end

function ossa_update!(c::Coupling, spcs::Vector{Species}, rxns::Vector{Reaction}, param, intensity::Float64, g::LightGraphs.DiGraph)
	τ = rand(Exponential(1 / intensity))
	jump = intensity * rand()
	μ = sample(rxns, jump)
	μ > 0 ? update!(spcs, rxns[μ]) : error("No reaction occurred!")

	intensity = update_propensities!(c, rxns, spcs, param, g)
	return τ, intensity
end

function update_propensities!(c::Loose, rxns, spcs, param)
	dependents = neighbors(g, μ)

	@inbounds for α in dependents
		intensity = intensity - rxns[α].propensity
		propensity!(rxns[α], spcs, param);
		intensity = intensity + rxns[α].propensity
	end
	return intensity
end

function update_propensities!(c::Tight, rxns, spcs, param, g)
	return compute_propensities!(rxns, spcs, param)
end

function presimulate!(spcs, rxns, params, init, n, itr)
	events = zeros(Float64, length(rxns))

	for i = 1:itr
		reset!(spcs, init)
		for k = 1:n
			intensity = compute_propensities!(rxns, spcs, params)
			τ = rand(Exponential(1 / intensity))
			u = rand()
		  jump = intensity * u
		  j = sample(rxns, jump)
		  j > 0 ? update!(spcs, rxns[j]) : error("No reaction occurred!")
			events[j] = events[j] + 1
		end
	end

	for i in eachindex(rxns)
		rxns[i].propensity = events[i]
	end

	return rxns
end

function init_ossa!(spcs, rxns, params, init, n, itr)
	rxns = presimulate!(spcs, rxns, params, init, n, itr)
	sort!(rxns, alg=Base.MergeSort, lt=isless, rev=true)
	return rxns
end

function isless(x::Reaction, y::Reaction)
	return x.propensity < y.propensity
end
