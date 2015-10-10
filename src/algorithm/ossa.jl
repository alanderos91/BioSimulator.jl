function ossa(model::Simulation, t_final::Float64, output::OutputType, dt::Float64, itr::Int)
	# Unpack model
	init = model.initial
	spcs = model.state
	rxns = model.rxns
	params = model.param

	n = round(Int, t_final / dt) + 1

	job = SimulationJob(itr)

	g = init_dep_graph(rxns)
	pq = init_apq(rxns)

	for i = 1:itr
		reset!(spcs, init)

		result = init_sr(output, spcs, n)
		ssa_steps = 0

		t = 0.0
		t_next = 0.0
		j = 1

		# Compute propensities and initialize the priority queue with propensities
		intensity = compute_propensities!(rxns, spcs, params)
		init_apq!(pq, rxns)

		while t < t_final
			t_next, j = update!(output, result, n, t, t_next, dt, j, spcs)
			τ, intensity = ossa_update!(spcs, rxns, params, intensity, g, pq)
			t = t + τ
			ssa_steps = ssa_steps + 1
		end
		update!(output, result, n, t_next, dt, j, spcs)
		job[i] = result
	end
	return job
end

function init_apq(rxns::Vector{Reaction})
	pq = Collections.PriorityQueue(Int,Float64,Base.Reverse)
	for j in eachindex(rxns)
		pq[j] = 0.0
	end
	return pq
end

function init_apq!(pq, rxns::Vector{Reaction})
	@inbounds for j in eachindex(rxns)
		pq[j] = rxns[j].propensity
	end
end

function ossa_update!(spcs::Vector{Species}, rxns::Vector{Reaction}, param, intensity::Float64, g::LightGraphs.DiGraph, pq)
	τ = rand(Exponential(1 / intensity))
	jump = intensity * rand()
	μ = sample(pq, jump)
	μ > 0 ? update!(spcs, rxns[μ]) : error("No reaction occurred!")

	# Update propensities
	dependents = neighbors(g, μ)

	@inbounds for α in dependents
		intensity = intensity - rxns[α].propensity
		propensity!(rxns[α], spcs, param); pq[α] = rxns[α].propensity
		intensity = intensity + rxns[α].propensity
	end
	return τ, intensity
end

function sample(pq, jump)
	ss = 0.0
	# This will iterate through the priority queue in order of decreasing
	# elements. The value i is the i-th largest element, j is the index
	# of the corresponding reaction, and a is the i-th largest propensity.
	for i = 1:length(pq)
		j, a = pq.xs[i]
		ss = ss + a
		if ss >= jump
			return j
		end
	end
	return 0
end
