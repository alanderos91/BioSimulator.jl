export simulate, Explicit, Uniform

abstract OutputType

immutable Explicit <: OutputType end

function update!(::Explicit, result, n, t, t_next, dt, j, spcs)
	update!(result, t, spcs)
	return t, j
end

function update!(::Explicit, result, n, t_next, dt, j, spcs)
	update!(result, t_next, spcs)
	return;
end

function init_sr(::Explicit, spcs, n)
	return SimulationResult(spcs)
end

immutable Uniform <: OutputType end

function update!(::Uniform, result, n, t, t_next, dt, j, spcs)
	while t >= t_next
		update!(result, t_next, spcs, j)
		j = j + 1
		t_next = t_next + dt
		if j > n; return t_next, j; end
	end
	return t_next, j
end

function update!(::Uniform, result, n, t_next, dt, j, spcs)
	while j <= n
		update!(result, t_next, spcs, j)
		j = j + 1
		t_next = t_next + 1
	end
	return;
end

function init_sr(::Uniform, spcs, n)
	return SimulationResult(spcs, n)
end

function simulate(model::Simulation, t_final::Float64, with::Symbol=:ssa; o::OutputType=Uniform(), dt::Float64=1.0, itr::Int=1, args...)
	d = Dict{Symbol, Float64}()
	for (key, value) in args
		d[key] = value
	end

	if with == :ssa
		ssa(model, t_final, o, dt, itr)
	elseif with == :ossa
		c = haskey(d, :c) ? d[:c] : Tight()
		steps = haskey(d, :steps) ? d[:steps] : 100
		samples = haskey(d, :samples) ? d[:samples] : 1
		ossa(model, t_final, o, dt, itr, c, steps, samples)
	elseif with == :frm
		frm(model, t_final, o, dt, itr)
	elseif with == :nrm
		nrm(model, t_final, o, dt, itr)
	elseif with == :sal
		tol = haskey(d, :tol) ? d[:tol] : 0.125
		thrsh = haskey(d, :thrsh) ? d[:thrsh] : 100.0
		ctrct = haskey(d, :ctrct) ? d[:ctrct] : 0.75
		sal(model, t_final, o, dt, itr, tol, thrsh, ctrct)
	end
end
