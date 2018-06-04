struct Direct end

function build_algorithm(::Direct, end_time, track_stats; kwargs...)
  SSA(end_time, track_stats)
end

struct FirstReaction end

function build_algorithm(::FirstReaction, end_time, track_stats; kwargs...)
  FRM(end_time, track_stats)
end

struct NextReaction end

function build_algorithm(::NextReaction, end_time, track_stats; kwargs...)
  NRM(end_time, track_stats)
end

struct OptimizedDirect end

function build_algorithm(::OptimizedDirect, end_time, track_stats; kwargs...)
  ODM(end_time, track_stats)
end

struct TauLeaping end

function build_algorithm(::TauLeaping, end_time, track_stats; ϵ::Float64=0.125, δ::Float64=2.0, β::Float64=0.75, kwargs...)
  OTL(end_time, ϵ, δ, β, track_stats)
end

struct StepAnticipation end

function build_algorithm(::StepAnticipation, end_time, track_stats; ϵ::Float64=0.125, δ::Float64=2.0, β::Float64=0.75, kwargs...)
  SAL(end_time, ϵ, δ, β, track_stats)
end