abstract Algorithm

immutable SSA <: Algorithm end
immutable FRM <: Algorithm end
immutable NRM <: Algorithm end
immutable SAL <: Algorithm end

abstract OutputType

immutable Explicit <: OutputType end
immutable Uniform <: OutputType end

function simulate(model::Simulation, t_final::Float64, with::Algorithm=SSA; output::OutputType=Uniform, dt::Float64=1.0, itr::Int=1)
	if with <: SSA
		ssa(model, t_final, output=output, dt=dt, itr=itr)
	elseif with <: FRM
		frm(model, t_final, output=output, dt=dt, itr=itr)
	elseif with <: NRM
		nrm(model, t_final, output=output, dt=dt, itr=itr)
	elseif with <: SAL
		sal(model, t_final, output=output, dt=dt, itr=itr)
	end
end

export simulate, SSA, FRM, NRM, SAL, Explicit, Uniform
