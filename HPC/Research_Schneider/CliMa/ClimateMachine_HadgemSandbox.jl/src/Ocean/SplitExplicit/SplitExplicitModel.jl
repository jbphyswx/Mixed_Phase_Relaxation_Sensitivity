module SplitExplicit

using StaticArrays

using ..HydrostaticBoussinesq
using ..ShallowWater

using ...VariableTemplates
using ...MPIStateArrays
using ...Mesh.Geometry
using ...DGMethods
using ...BalanceLaws

include("VerticalIntegralModel.jl")

import ...BalanceLaws:
    initialize_states!,
    tendency_from_slow_to_fast!,
    cummulate_fast_solution!,
    reconcile_from_fast_to_slow!

@inline initialize_states!(
    slow::HydrostaticBoussinesqModel,
    fast::ShallowWaterModel,
    _...,
) = nothing

@inline tendency_from_slow_to_fast!(
    slow::HydrostaticBoussinesqModel,
    fast::ShallowWaterModel,
    _...,
) = nothing

@inline cummulate_fast_solution!(
    slow::HydrostaticBoussinesqModel,
    fast::ShallowWaterModel,
    _...,
) = nothing

@inline reconcile_from_fast_to_slow!(
    slow::HydrostaticBoussinesqModel,
    fast::ShallowWaterModel,
    _...,
) = nothing

end
