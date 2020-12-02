##### Moisture tendencies

#####
##### First order fluxes
#####

function flux(::Advect{TotalMoisture}, m, state, aux, t, ts, direction)
    u = state.ρu / state.ρ
    return u * state.moisture.ρq_tot
end

function flux(::Advect{LiquidMoisture}, m, state, aux, t, ts, direction)
    u = state.ρu / state.ρ
    return u * state.moisture.ρq_liq
end

function flux(::Advect{IceMoisture}, m, state, aux, t, ts, direction)
    u = state.ρu / state.ρ
    return u * state.moisture.ρq_ice
end

#####
##### Sources
#####

function source(
    s::Subsidence{TotalMoisture},
    m,
    state,
    aux,
    t,
    ts,
    direction,
    diffusive,
)
    z = altitude(m, aux)
    w_sub = subsidence_velocity(s, z)
    k̂ = vertical_unit_vector(m, aux)
    return -state.ρ * w_sub * dot(k̂, diffusive.moisture.∇q_tot)
end

export CreateClouds
"""
    CreateClouds{PV <: Union{LiquidMoisture,IceMoisture}} <: TendencyDef{Source, PV}

A source/sink to `q_liq` and `q_ice` implemented as a relaxation towards
equilibrium in the Microphysics module.
The default relaxation timescales are defined in CLIMAParameters.jl.
"""
struct CreateClouds{PV <: Union{LiquidMoisture, IceMoisture}} <:
       TendencyDef{Source, PV} end

CreateClouds() = (CreateClouds{LiquidMoisture}(), CreateClouds{IceMoisture}())

function source(
    s::CreateClouds{LiquidMoisture},
    m,
    state,
    aux,
    t,
    ts,
    direction,
    diffusive,
)
    # get current temperature and phase partition
    FT = eltype(state)
    q = PhasePartition(ts)
    q_sat_liq = q_vap_saturation_liquid(ts)
    q_vap = vapor_specific_humidity(q)

    # cloud condensate as relaxation source terms
    S_q_liq = conv_q_vap_to_q_liq_ice2(m.param_set.microphys.liq, q_sat_liq, q_vap)

    return state.ρ * S_q_liq
end

function source(
    s::CreateClouds{IceMoisture},
    m,
    state,
    aux,
    t,
    ts,
    direction,
    diffusive,
)
    FT = eltype(state)
    q = PhasePartition(ts)
    q_sat_ice = q_vap_saturation_ice(ts)
    q_vap = vapor_specific_humidity(q)

    # cloud condensate as relaxation source terms
    S_q_ice = conv_q_vap_to_q_liq_ice2(m.param_set.microphys.ice, q_sat_ice, q_vap)

    return state.ρ * S_q_ice
end

function source(
    s::RemovePrecipitation{TotalMoisture},
    m,
    state,
    aux,
    t,
    ts,
    direction,
    diffusive,
)
    if has_condensate(ts)
        nt = compute_precip_params(s, aux, ts)
        return state.ρ * nt.S_qt
    else
        FT = eltype(state)
        return FT(0)
    end
end

function source(
    s::Rain_1M{TotalMoisture},
    m,
    state,
    aux,
    t,
    ts,
    direction,
    diffusive,
)
    nt = compute_rain_params(m, state, aux, t, ts)
    return state.ρ * nt.S_qt
end

function source(
    s::Rain_1M{LiquidMoisture},
    m,
    state,
    aux,
    t,
    ts,
    direction,
    diffusive,
)
    nt = compute_rain_params(m, state, aux, t, ts)
    return state.ρ * nt.S_ql
end
