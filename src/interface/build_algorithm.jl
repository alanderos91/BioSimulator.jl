struct Direct end

function build_algorithm(::Direct, end_time; kwargs...)
  SSA(end_time)
end

struct FirstReaction end

function build_algorithm(::FirstReaction, end_time; kwargs...)
  FRM(end_time)
end

struct NextReaction end

function build_algorithm(::NextReaction, end_time; kwargs...)
  NRM(end_time)
end

struct OptimizedDirect end

function build_algorithm(::OptimizedDirect, end_time; kwargs...)
  ODM(end_time)
end

struct TauLeaping end

function build_algorithm(::TauLeaping, end_time; ϵ::Float64=0.125, δ::Float64=100.0, β::Float64=0.75, kwargs...)
  OTL(end_time, ϵ, δ, β)
end

struct StepAnticipation end

function build_algorithm(::StepAnticipation, end_time; ϵ::Float64=0.125, δ::Float64=100.0, β::Float64=0.75, kwargs...)
  SAL(end_time, ϵ, δ, β)
end