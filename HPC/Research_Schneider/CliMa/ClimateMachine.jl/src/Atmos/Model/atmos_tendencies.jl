#####
##### Tendency specification
#####

import ..BalanceLaws: eq_tends

#####
##### Sources
#####

# --------- Some of these methods are generic or
#           temporary during transition to new specification:
filter_source(pv::PrognosticVariable, s) = nothing
# Sources that have been added to new specification:
filter_source(pv::PV, s::Subsidence{PV}) where {PV <: PrognosticVariable} = s
filter_source(pv::PV, s::Gravity{PV}) where {PV <: Momentum} = s
filter_source(pv::PV, s::GeostrophicForcing{PV}) where {PV <: Momentum} = s
filter_source(pv::PV, s::Coriolis{PV}) where {PV <: Momentum} = s
filter_source(pv::PV, s::RayleighSponge{PV}) where {PV <: Momentum} = s
filter_source(pv::PV, s::CreateClouds{PV}) where {PV <: LiquidMoisture} = s
filter_source(pv::PV, s::CreateClouds{PV}) where {PV <: IceMoisture} = s

# Filter sources / empty elements
filter_sources(t::Tuple) = filter(x -> !(x == nothing), t)
filter_sources(pv::PrognosticVariable, srcs) =
    filter_sources(map(s -> filter_source(pv, s), srcs))

# Entry point
eq_tends(pv::PrognosticVariable, m::AtmosModel, ::Source) =
    filter_sources(pv, m.source)
# ---------

#####
##### First order fluxes
#####

# Mass
eq_tends(pv::PV, ::AtmosModel, ::Flux{FirstOrder}) where {PV <: Mass} =
    (Advect{PV}(),)

# Momentum
eq_tends(pv::PV, m::AtmosModel, ::Flux{FirstOrder}) where {PV <: Momentum} =
    (Advect{PV}(), PressureGradient{PV}())

# Energy
eq_tends(pv::PV, ::AtmosModel, ::Flux{FirstOrder}) where {PV <: Energy} =
    (Advect{PV}(), Pressure{PV}())

# Moisture
eq_tends(pv::PV, ::AtmosModel, ::Flux{FirstOrder}) where {PV <: Moisture} =
    (Advect{PV}(),)

#####
##### Second order fluxes
#####

# Mass
moist_diffusion(pv::PV, ::DryModel) where {PV <: Mass} = ()
moist_diffusion(pv::PV, ::MoistureModel) where {PV <: Mass} =
    (MoistureDiffusion{PV}(),)
eq_tends(pv::PV, m::AtmosModel, ::Flux{SecondOrder}) where {PV <: Mass} =
    (moist_diffusion(pv, m.moisture)...,)

# Momentum
eq_tends(pv::PV, ::AtmosModel, ::Flux{SecondOrder}) where {PV <: Momentum} =
    (ViscousStress{PV}(),)

# Energy
eq_tends(pv::PV, ::AtmosModel, ::Flux{SecondOrder}) where {PV <: Energy} =
    (ViscousProduction{PV}(), EnthalpyProduction{PV}())

# Moisture
eq_tends(pv::PV, ::AtmosModel, ::Flux{SecondOrder}) where {PV <: Moisture} = ()
