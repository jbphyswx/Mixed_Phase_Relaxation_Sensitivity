export GCMModel, NoGCM, HadGEMVertical, CMIP_cfsite_Vertical

abstract type GCMModel end

vars_state(::GCMModel, ::Prognostic, FT) = @vars()
vars_state(::GCMModel, ::Auxiliary, FT) = @vars()
vars_state(::GCMModel, ::Gradient, FT) = @vars()
vars_state(::GCMModel, ::GradientFlux, FT) = @vars()
vars_state(::GCMModel, ::UpwardIntegrals, FT) = @vars()
vars_state(::GCMModel, ::DownwardIntegrals, FT) = @vars()

function atmos_nodal_update_auxiliary_state!(
    ::GCMModel,
    ::AtmosModel,
    state::Vars,
    aux::Vars,
    t::Real,
) end
function integral_set_auxiliary_state!(
    ::GCMModel,
    integ::Vars,
    state::Vars,
    aux::Vars,
) end
function integral_load_auxiliary_state!(::GCMModel, aux::Vars, integ::Vars) end
function reverse_integral_set_auxiliary_state!(
    ::GCMModel,
    integ::Vars,
    state::Vars,
    aux::Vars,
) end
function reverse_integral_load_auxiliary_state!(
    ::GCMModel,
    aux::Vars,
    integ::Vars,
) end
function flux_radiation!(
    ::GCMModel,
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) end
@inline function atmos_nodal_update_auxiliary_state!(
    gcminfo::GCMModel,
    atmos::AtmosModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    nothing
end

function compute_gradient_argument!(
    gcminfo::GCMModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    nothing
end

function compute_gradient_flux!(
    gcminfo::GCMModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    nothing
end

function flux_first_order_gcm!(
    gcminfo::GCMModel,
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    nothing
end

function flux_second_order!(
    gcminfo::GCMModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    D_t,
)
    nothing
end
function flux_second_order!(gcminfo::GCMModel, flux::Grad, state::Vars, d_q_tot)
    nothing
end

struct NoGCM <: GCMModel end

"""
    Container for GCM variables from HadGEM2-A forcing,
used in the AMIP experiments
"""
struct HadGEMVertical <: GCMModel end

vars_state(m::HadGEMVertical, ::Auxiliary, FT) = @vars(
    ta::FT,
    hus::FT,
    ua::FT,
    va::FT,
    tntΣhava::FT,
    Σtemp_tendency::FT,
    Σqt_tendency::FT,
    w_s::FT,
)

vars_state(::HadGEMVertical, ::Prognostic, FT) = @vars()
vars_state(::HadGEMVertical, ::Gradient, FT) = @vars(ta::FT, hus::FT)
vars_state(::HadGEMVertical, ::GradientFlux, FT) = @vars(∇ᵥta::FT, ∇ᵥhus::FT)

@inline function atmos_nodal_update_auxiliary_state!(
    gcminfo::HadGEMVertical,
    atmos::AtmosModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    nothing
end

function compute_gradient_argument!(
    gcminfo::HadGEMVertical,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    transform.gcminfo.ta = aux.gcminfo.ta
    transform.gcminfo.hus = aux.gcminfo.hus
end

function compute_gradient_flux!(
    gcminfo::HadGEMVertical,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    diffusive.gcminfo.∇ᵥta = ∇transform.gcminfo.ta[3]
    diffusive.gcminfo.∇ᵥhus = ∇transform.gcminfo.hus[3]
end

function flux_first_order_gcm!(
    gcminfo::HadGEMVertical,
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    nothing
end

function flux_second_order!(
    gcminfo::HadGEMVertical,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    D_t,
)
    nothing
end
function flux_second_order!(
    gcminfo::HadGEMVertical,
    flux::Grad,
    state::Vars,
    d_q_tot,
)
    nothing
end













"""
    Container for GCM variables from models for forcing my cfsite runs,
used in the AMIP(0-4K) experiments
"""
struct CMIP_cfsite_Vertical <: GCMModel end

vars_state(m::CMIP_cfsite_Vertical, ::Auxiliary, FT) = @vars(
    ta::FT,
    hus::FT,
    ua::FT,
    va::FT,
    tntΣhava::FT,
    Σtemp_tendency::FT,
    Σqt_tendency::FT,
    w_s::FT,
    ρ::FT,
    ρe::FT,
    ρq_tot::FT,
    ρq_rai::FT,
    ρq_sno::FT, # ion think this one is used yet but future prufe my g

)



vars_state(::CMIP_cfsite_Vertical, ::Prognostic, FT) = @vars()
vars_state(::CMIP_cfsite_Vertical, ::Gradient, FT) = @vars(ta::FT, hus::FT)
vars_state(::CMIP_cfsite_Vertical, ::GradientFlux, FT) =
    @vars(∇ᵥta::FT, ∇ᵥhus::FT)

@inline function atmos_nodal_update_auxiliary_state!(
    gcminfo::CMIP_cfsite_Vertical,
    atmos::AtmosModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    nothing
end

function compute_gradient_argument!(
    gcminfo::CMIP_cfsite_Vertical,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    transform.gcminfo.ta = aux.gcminfo.ta
    transform.gcminfo.hus = aux.gcminfo.hus
end

function compute_gradient_flux!(
    gcminfo::CMIP_cfsite_Vertical,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    diffusive.gcminfo.∇ᵥta = ∇transform.gcminfo.ta[3]
    diffusive.gcminfo.∇ᵥhus = ∇transform.gcminfo.hus[3]
end

function flux_first_order_gcm!(
    gcminfo::CMIP_cfsite_Vertical,
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    nothing
end

function flux_second_order!(
    gcminfo::CMIP_cfsite_Vertical,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    D_t,
)
    nothing
end
function flux_second_order!(gcminfo::CMIP_cfsite_Vertical, flux::Grad, state::Vars, d_q_tot)
    nothing
end