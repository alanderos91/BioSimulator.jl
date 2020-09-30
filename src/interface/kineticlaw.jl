##### general kinetic laws #####

abstract type KineticLaw end
struct MassAction  <: KineticLaw end
struct DefaultIPSLaw <: KineticLaw end