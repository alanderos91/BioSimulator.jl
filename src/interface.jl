export simulate, Explicit, Uniform

abstract OutputType

immutable Explicit <: OutputType end
immutable Uniform <: OutputType end

function update!(::Explicit, result, n, t, t_next, dt, j, spcs)
	update!(result, t, spcs)
	return t_next, j
end

function init_sr(::Explicit, spcs, n)
	return SimulationResult(spcs)
end

function update!(::Uniform, result, n, t, t_next, dt, j, spcs)
	while t >= t_next
		update!(result, t_next, spcs, j)
		j = j + 1
		t_next = t_next + dt
		if j > n; return t_next, j; end
	end
	return t_next, j
end

function init_sr(::Uniform, spcs, n)
	return SimulationResult(spcs, n)
end

function simulate(model::Simulation, t_final::Float64, with::Symbol=:ssa; o::OutputType=Uniform(), dt::Float64=1.0, itr::Int=1)
	if with == :ssa
		ssa(model, t_final, o, dt, itr)
	elseif with == :frm
		frm(model, t_final, o, dt, itr)
	elseif with == :nrm
		nrm(model, t_final, o, dt, itr)
	elseif with == :sal
		sal(model, t_final, o, dt, itr)
	end
end
