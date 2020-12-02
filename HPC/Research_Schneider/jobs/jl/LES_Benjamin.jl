#!/usr/bin/env julia --project

# # Southern Ocean Mixed Phase cfSites simulations 

# Adapted from:
# BOMEX cases in Gulf - https://github.com/CliMA/ClimateMachine.jl/blob/master/experiments/AtmosLES/bomex_les.jl)
# Zhaoyi branch       - https://github.com/CliMA/ClimateMachine.jl/tree/as/hadgem-sandbox

# General Julia Modules
using ArgParse
using Distributions
using DocStringExtensions
using LinearAlgebra
using Printf
using Random
using StaticArrays
using Test
using Dierckx
using NCDatasets

# Climate Machine Modules
using ClimateMachine
using ClimateMachine.Atmos
using ClimateMachine.Orientations
using ClimateMachine.ConfigTypes
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.Diagnostics
using ClimateMachine.GenericCallbacks
using ClimateMachine.Mesh.Filters
using ClimateMachine.Mesh.Grids
using ClimateMachine.ODESolvers
using ClimateMachine.Thermodynamics
using ClimateMachine.TurbulenceClosures
using ClimateMachine.TurbulenceConvection
using ClimateMachine.VariableTemplates
using ClimateMachine.BalanceLaws:
    BalanceLaw, Auxiliary, Gradient, GradientFlux, Prognostic

using CLIMAParameters
using CLIMAParameters.Planet: e_int_v0, grav, day
using CLIMAParameters.Atmos.Microphysics
using ClimateMachine.Atmos: altitude, recover_thermo_state


# ---- Physics specific imports 
import ClimateMachine.DGMethods: vars_state_conservative, vars_state_auxiliary
import ClimateMachine.Atmos: source!, atmos_source!, altitude
import ClimateMachine.Atmos: compute_gradient_flux!, thermo_state

struct LiquidParameterSet <: AbstractLiquidParameterSet end
struct IceParameterSet <: AbstractIceParameterSet end

struct MicropysicsParameterSet{L, I} <: AbstractMicrophysicsParameterSet
    liq::L
    ice::I
end

struct EarthParameterSet{M} <: AbstractEarthParameterSet
    microphys::M
end
microphys = MicropysicsParameterSet(LiquidParameterSet(), IceParameterSet())
const param_set = EarthParameterSet(microphys)


struct GCMRelaxation{FT} <: Source
    τ_relax::FT
end
function atmos_source!(
    s::GCMRelaxation,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    return nothing
end

# Temperature tendency term, ∂T∂t (from hadgem branch?)
"""
    TemperatureTendency <: Source
Temperature tendency for the LES configuration based on quantities 
from a GCM. Quantities are included in standard CMIP naming format. 
Tendencies included here are 
    tntha = temperature tendency due to horizontal advection
    tntva = temperature tendency due to vertical advection
    tntr = temperature tendency due to radiation fluxes
    ∂T∂z = temperature vertical gradient from GCM values
"""
struct TemperatureTendency <: Source end
function atmos_source!(
    s::TemperatureTendency,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction
)
    # Establish problem float-type
    FT = eltype(state)
    _grav = grav(atmos.param_set)
    # Establish vertical orientation
    k̂ = vertical_unit_vector(atmos, aux)
    _e_int_v0 = e_int_v0(atmos.param_set)
    # Unpack vertical gradients
    ∂qt∂z = dot(diffusive.moisture.∇q_tot_gcm, k̂)
    ∂T∂z = dot(diffusive.moisture.∇T_gcm, k̂)
    w_adv = -aux.gcminfo.wap / (aux.gcminfo.ρ * _grav)
    # Establish thermodynamic state
    TS = thermo_state(atmos, state, aux)
    cvm = cv_m(TS)
    # Compute tendency terms
    # -- Temperature contribution
    source.ρe += cvm * state.ρ * aux.gcminfo.tntha
    source.ρe += cvm * state.ρ * aux.gcminfo.tntva
    source.ρe += cvm * state.ρ * aux.gcminfo.tntr
    source.ρe += cvm * state.ρ * ∂T∂z * w_adv
    # -- Moisture contribution
    source.ρe += _e_int_v0 * state.ρ * aux.gcminfo.tnhusha
    source.ρe += _e_int_v0 * state.ρ * aux.gcminfo.tnhusva
    source.ρe += _e_int_v0 * state.ρ * ∂qt∂z * w_adv
    # GPU-friendly return nothing
    return nothing
end


# Moisture tendency ∂qt∂t
"""
    MoistureTendency <: Source 
Moisture tendency for the LES configuration based on quantities 
from a GCM. Quantities are included in standard CMIP naming format. 
Tendencies included here are 
    tnhusha = moisture tendency due to horizontal advection
    tnhusva = moisture tendency due to vertical advection
    ∂qt∂z = moisture vertical gradient from GCM values
"""
struct MoistureTendency <: Source end
function atmos_source!(
    s::MoistureTendency,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction
)
    # Establish problem float-type
    FT = eltype(state)
    _grav = grav(atmos.param_set)
    k̂ = vertical_unit_vector(atmos, aux)
    # Establish vertical orientation
    ∂qt∂z = dot(diffusive.moisture.∇q_tot_gcm, k̂)
    w_adv = -aux.gcminfo.wap / (aux.gcminfo.ρ * _grav)
    # Establish thermodynamic state
    TS = thermo_state(atmos, state, aux)
    cvm = cv_m(TS)
    # Compute tendency terms
    source.moisture.ρq_tot += state.ρ * aux.gcminfo.tnhusha
    source.moisture.ρq_tot += state.ρ * aux.gcminfo.tnhusva
    source.moisture.ρq_tot += state.ρ * ∂qt∂z * w_adv
    
    source.ρ += state.ρ * aux.gcminfo.tnhusha
    source.ρ += state.ρ * aux.gcminfo.tnhusva
    source.ρ += state.ρ * ∂qt∂z * w_adv
    # GPU-friendly return nothing
    return nothing
end



# Large-scale subsidence forcing
"""
    Subsidence <: Source 
Subsidence tendency, given a vertical velocity at the large scale, 
obtained from the GCM data. 
    wap = GCM vertical velocity [Pa s⁻¹]. Note the conversion required
"""
struct SubsidenceTendency <: Source end
function atmos_source!(
    s::SubsidenceTendency,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction
)
    _grav = grav(atmos.param_set)
    # Establish vertical orientation
    k̂ = vertical_unit_vector(atmos, aux)
    # Establish subsidence velocity
    w_s = -aux.gcminfo.wap / aux.gcminfo.ρ / _grav
    # Compute tendency terms
    source.ρe -= state.ρ * w_s * dot(k̂, diffusive.∇h_tot)
    source.moisture.ρq_tot -= state.ρ * w_s * dot(k̂, diffusive.moisture.∇q_tot)
    # GPU-friendly return nothing
    return nothing
end



# Sponge relaxation
"""
  LinearSponge (Source)
Two parameter sponge (α_max, γ) for velocity relaxation to a reference
state.
    α_max = Sponge strength (can be interpreted as timescale)
    γ = Sponge exponent
    z_max = Domain height
    z_sponge = Altitude at which sponge layer starts
"""
struct LinearSponge{FT} <: Source
    "Maximum domain altitude (m)"
    z_max::FT
    "Altitude at with sponge starts (m)"
    z_sponge::FT
    "Sponge Strength 0 ⩽ α_max ⩽ 1"
    α_max::FT
    "Sponge exponent"
    γ::FT
end
function atmos_source!(
    s::LinearSponge,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction
)
    #Unpack sponge parameters
    FT = eltype(state)
    z_max = s.z_max
    z_sponge = s.z_sponge
    α_max = s.α_max
    γ = s.γ
    # Establish sponge relaxation velocity
    u_geo = SVector(aux.gcminfo.ua, aux.gcminfo.va, 0)
    z = altitude(atmos, aux)
    # Accumulate sources
    if z_sponge <= z
        r = (z - z_sponge) / (z_max - z_sponge)
        β_sponge = α_max .* sinpi(r/FT(2)) * sinpi(r/FT(2)) * sinpi(r/ FT(2)) * sinpi(r/FT(2))#.^ γ
        source.ρu -= β_sponge * (state.ρu .- state.ρ * u_geo)
    end
    # GPU-friendly return nothing
    return nothing
end


"""


Insert something here to base initialization off of an input file...
Akshay and Zhaoyi had a custom file 


"""

# Initialise the CFSite experiment :D! 
const seed = MersenneTwister(0) #??? what does this do
function init_cfsites!(problem, bl, state, aux, (x, y, z), t, spl)
    FT = eltype(state)
    _grav = grav(bl.param_set)

    # Unpack splines, interpolate to z coordinate at 
    # present grid index. (Functions are all pointwise)
    ta = FT(spl.spl_ta(z))
    hus = FT(spl.spl_hus(z))
    ua = FT(spl.spl_ua(z))
    va = FT(spl.spl_va(z))
    pfull = FT(spl.spl_pfull(z))
    tntha = FT(spl.spl_tntha(z))
    tntva = FT(spl.spl_tntva(z))
    tntr = FT(spl.spl_tntr(z))
    tnhusha = FT(spl.spl_tnhusha(z))
    tnhusva = FT(spl.spl_tnhusva(z))
    wap = FT(spl.spl_wap(z))
    ρ_gcm = FT(1 / spl.spl_alpha(z))

    # Compute field properties based on interpolated data
    ρ = air_density(bl.param_set, ta, pfull, PhasePartition(hus))
    e_int = internal_energy(bl.param_set, ta, PhasePartition(hus))
    e_kin = (ua^2 + va^2) / 2
    e_pot = _grav * z
    # Assignment of state variables
    state.ρ = ρ
    state.ρu = ρ * SVector(ua, va, 0)
    state.ρe = ρ * (e_kin + e_pot + e_int)
    state.moisture.ρq_tot = ρ * hus
    if z <= FT(400)
        state.ρe += rand(seed) * FT(1 / 100) * (state.ρe)
        state.moisture.ρq_tot +=
            rand(seed) * FT(1 / 100) * (state.moisture.ρq_tot)
    end

    # Assign and store the ref variable for sources
    aux.gcminfo.ρ = ρ_gcm
    aux.gcminfo.pfull = pfull
    aux.gcminfo.ta = ta
    aux.gcminfo.hus = hus
    aux.gcminfo.tntha = tntha
    aux.gcminfo.tntva = tntva
    aux.gcminfo.ua = ua
    aux.gcminfo.va = va
    aux.gcminfo.tntr = tntr
    aux.gcminfo.tnhusha = tnhusha
    aux.gcminfo.tnhusva = tnhusva
    aux.gcminfo.wap = wap

    return nothing
end




function config_cfsites(
    FT,
    N,
    resolution,
    xmax,
    ymax,
    zmax,
    hfls,
    hfss,
    ts,
    groupid,
)
    # Boundary Conditions
    u_star = FT(0.28)

    problem = AtmosProblem(
        boundarycondition = (
            AtmosBC(
                momentum = Impenetrable(DragLaw(
                    (state, aux, t, normPu_int) -> (u_star / normPu_int)^2,
                )),
                energy = PrescribedEnergyFlux((state, aux, t) -> hfls + hfss),
                moisture = PrescribedMoistureFlux(
                    (state, aux, t) ->
                        hfls / latent_heat_vapor(param_set, ts),
                ),
            ),
            AtmosBC(),
        ),
        init_state_prognostic = init_cfsites!,
    )
    model = AtmosModel{FT}(
        AtmosLESConfigType,
        param_set;
        problem = problem,
        turbulence = Vreman{FT}(0.23),
        source = (
            Gravity(),
            LinearSponge{FT}(zmax, zmax * 0.85, 1, 4),
            MoistureTendency(),
            EnergyTendency(),
            SubsidenceTendency(),
        ),
        moisture = EquilMoist{FT}(; maxiter = 5, tolerance = FT(2)),
        #ZS: hyperdiffusion?
        #hyperdiffusion = DryBiharmonic{FT}(12*3600),
        gcminfo = HadGEM(),
    )

    # Timestepper options

    # Explicit Solver
    ex_solver = ClimateMachine.ExplicitSolverType()

    # Multirate Explicit Solver
    mrrk_solver = ClimateMachine.MultirateSolverType(
        fast_model = AtmosAcousticGravityLinearModel,
        slow_method = LSRK144NiegemannDiehlBusch,
        fast_method = LSRK144NiegemannDiehlBusch,
        timestep_ratio = 10,
    )

    # IMEX Solver Type
    imex_solver = ClimateMachine.IMEXSolverType()

    # Configuration
    config = ClimateMachine.AtmosLESConfiguration(
        forcingfile * "_$groupid",
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        param_set,
        init_cfsites!,
        #ZS: multi-rate?
        solver_type = imex_solver,
        model = model,
    )
    return config
end

# Define the diagnostics configuration (Atmos-Default)
function config_diagnostics(driver_config)
    default_dgngrp = setup_atmos_default_diagnostics(
        AtmosLESConfigType(),
        "2500steps",
        driver_config.name,
    )
    core_dgngrp = setup_atmos_core_diagnostics(
        AtmosLESConfigType(),
        "2500steps",
        driver_config.name,
    )
    return ClimateMachine.DiagnosticsConfiguration([
        default_dgngrp,
        core_dgngrp,
    ])
end

function main()

    # Provision for custom command line arguments
    # Convenience args for slurm-array launches
    cfsite_args = ArgParseSettings(autofix_names = true)
    add_arg_group!(cfsite_args, "HadGEM2-A_SiteInfo")
    @add_arg_table! cfsite_args begin
        "--group-id"
        help = "Specify CFSite-ID for GCM data"
        metavar = "site<number>"
        arg_type = String
        default = "site17"
    end
    cl_args =
        ClimateMachine.init(parse_clargs = true, custom_clargs = cfsite_args)
    groupid = cl_args["group_id"]

    # Working precision
    FT = Float64
    # DG polynomial order
    N = 4
    # Domain resolution and size
    Δh = FT(75)
    Δv = FT(20)
    resolution = (Δh, Δh, Δv)
    # Domain extents
    xmax = FT(1800)
    ymax = FT(1800)
    zmax = FT(4000)
    # Simulation time
    t0 = FT(0)
    timeend = FT(600)
    #timeend = FT(3600 * 6)
    # Courant number
    CFL = FT(0.8)

    # Execute the get_gcm_info function
    (
        z,
        ta,
        hus,
        ua,
        va,
        pfull,
        tntha,
        tntva,
        tntr,
        tnhusha,
        tnhusva,
        wap,
        hfls,
        hfss,
        ts,
        alpha,
    ) = get_gcm_info(groupid)

    # Drop dimensions for compatibility with Dierckx
    z = dropdims(z; dims = 2)
    # Create spline objects and pass them into a named tuple
    splines = (
        spl_ta = Spline1D(z, view(ta, :, 1)),
        spl_pfull = Spline1D(z, view(pfull, :, 1)),
        spl_ua = Spline1D(z, view(ua, :, 1)),
        spl_va = Spline1D(z, view(va, :, 1)),
        spl_hus = Spline1D(z, view(hus, :, 1)),
        spl_tntha = Spline1D(z, view(tntha, :, 1)),
        spl_tntva = Spline1D(z, view(tntva, :, 1)),
        spl_tntr = Spline1D(z, view(tntr, :, 1)),
        spl_tnhusha = Spline1D(z, view(tnhusha, :, 1)),
        spl_tnhusva = Spline1D(z, view(tnhusva, :, 1)),
        spl_wap = Spline1D(z, view(wap, :, 1)),
        spl_alpha = Spline1D(z, view(alpha, :, 1)),
    )

    # Set up driver configuration
    driver_config = config_cfsites(
        FT,
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        hfls,
        hfss,
        ts,
        groupid,
    )
    # Set up solver configuration
    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        splines;
        init_on_cpu = true,
        Courant_number = CFL,
    )
    # Set up diagnostic configuration
    dgn_config = config_diagnostics(driver_config)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            ("moisture.ρq_tot",),
            solver_config.dg.grid,
            TMARFilter(),
        )
        nothing
    end

    #ZS: cutoff filter?
    filterorder = 2 * N
    filter = BoydVandevenFilter(solver_config.dg.grid, 0, filterorder)
    cbfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            AtmosFilterPerturbations(driver_config.bl),
            solver_config.dg.grid,
            filter,
        )
        nothing
    end

    cutoff_filter = CutoffFilter(solver_config.dg.grid, N - 1)
    cbcutoff = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(solver_config.Q, 1:6, solver_config.dg.grid, cutoff_filter)
        nothing
    end

    # Invoke solver (calls solve! function for time-integrator)
    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        #ZS: only tmar?
        user_callbacks = (cbtmarfilter,),
        check_euclidean_distance = true,
    )

end
main()

























