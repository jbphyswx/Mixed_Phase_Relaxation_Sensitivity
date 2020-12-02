""" 
-- TO DO -- 
> SWITCH TO ONLY ONE GCM PARSER, OUTSIDE THE MAIN() FCN AND THEN JUST CALL IT INSIDE MAIN TO AVOID HAVING TO DUPLICATE THE ENTRIRE CL_ARGS PARSING err_code
> CHECK HOW IMEX VS EX SOLVERS DIFFER AND WHAT MIGHT BE BEST TO USE...
> CHECK HOW THE MOISTURE RELAXATION WORKS TO BE SURE IT'S NOT ALSO THE LIQUID AND ICE THAT ARE BEING RELAXED
> PASS THROUGH SOLVER CHOICE TO OUTSIDE...

"""

# General Julia modules
using ArgParse
using Dierckx
using Distributions
using DocStringExtensions
using LinearAlgebra
using NCDatasets
using Printf
using Random
using StaticArrays
using Test

# using Pkg
# Pkg.add(CFTime)
# using CFTime
using  Dates # added by me for forcing fcn
import Statistics

# ClimateMachine Modules
using ClimateMachine

using ClimateMachine.Atmos
using ClimateMachine.ConfigTypes
using ClimateMachine.GenericCallbacks
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.Diagnostics
using ClimateMachine.Mesh.Filters
using ClimateMachine.Mesh.Grids
using ClimateMachine.Orientations
using ClimateMachine.Thermodynamics
using ClimateMachine.TurbulenceClosures
using ClimateMachine.TurbulenceConvection
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.Writers

using CLIMAParameters
using CLIMAParameters.Planet: R_d, planet_radius, grav, MSLP, molmass_ratio, e_int_v0, day, LH_v0
using CLIMAParameters.Atmos.Microphysics

# Physics specific imports ( is this still encessary?)
using ClimateMachine.Atmos: altitude, recover_thermo_state
import ClimateMachine.Atmos: source!, atmos_source!

# struct LiquidParameterSet <: AbstractLiquidParameterSet
    # τ_cond_evap:
# end
# struct IceParameterSet <: AbstractIceParameterSet
#     τ_sub_dep:: ice_param_set
# end

struct LiquidParameterSet <: AbstractLiquidParameterSet end
struct IceParameterSet    <: AbstractIceParameterSet    end
struct RainParameterSet   <: AbstractRainParameterSet   end
struct SnowParameterSet   <: AbstractSnowParameterSet   end

# struct MicrophysicsParameterSet{L, I} <: AbstractMicrophysicsParameterSet
#     liq::L
#     ice::I
# end
# const microphys = MicrophysicsParameterSet(LiquidParameterSet(), IceParameterSet()) # This const declaration is non-essential but hey why not

struct MicrophysicsParameterSet{L, I, R, S} <: AbstractMicrophysicsParameterSet # see dycoms example, https://github.com/CliMA/ClimateMachine.jl/blob/3cd5f471bbee32e8cee037e70bc36c0f35b05a5c/experiments/AtmosLES/dycoms.jl
    liq::L
    ice::I
    rai::R
    sno::S
end
const microphys = MicrophysicsParameterSet(LiquidParameterSet(), IceParameterSet(), RainParameterSet(), SnowParameterSet())

struct EarthParameterSet{M} <: AbstractEarthParameterSet
    microphys::M
end

const param_set = EarthParameterSet(microphys) # is this bad to do? It's being modified, no? 

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #


function parse_commandline(;parse_clargs=true, init_driver=false)

    # Working precision
    FT = Float64

    # Provision for custom command line arguments
    # Convenience args for slurm-array launches
    cfsite_args = ArgParseSettings(autofix_names = true)
    add_arg_group!(cfsite_args, "HadGEM2-A_SiteInfo")
    @add_arg_table! cfsite_args begin
        "--data_path"
            help = "Specify data path for loading model data from"
            arg_type = String
            # default = "/home/jbenjami/Research_Schneider/Data/cfsites/CMIP5/CFMIP2/forcing/"
        "--model"
            help = "Specify model data we're loading"
            arg_type = String
            # default = "HadGEM2-A" # or just nothing perhaps is safer?
        "--exper"
            help = "Specify the experiment we're loading from"
            arg_type = String
            # default = "amip"
        "--rip"
            help = "Specify the RIP (ensemble) value we're loading from"
            arg_type = String
            # default = "r1i1p1"
        "--years"
            help = "Specify the years we want to use data from"
            arg_type = String
            default = "all"
        "--months"
            help = "Specify the months we wish to use"
            arg_type = String
            default = "all"
        "--days"
            help = "Specify the days we wish to use"
            arg_type = String
            default = "all"
        "--hours"
            help = "Specify the hours we wish to use"
            arg_type = String
            default = "all"
        "--minutes"
            help = "Specify the minutes we wish to use"
            arg_type = String
            default = "all"
        "--seconds"
            help = "Specify the seconds we wish to use"
            arg_type = String
            default = "all"
        "--sites"
            # metavar = "site<number>"
            help = "Specify CFSite-IDs for GCM data (averages over selected sites)"
            arg_type = String
            # default = "all"
        "--delta_h"
            help = "Specify horizontal resolution (m)"
            arg_type = FT
            default = FT(75)
        "--delta_v"
            help = "Specify vertical resolution (m)"
            arg_type = FT
            default = FT(20)
        "--xmax"
            help = "Specify maximum x extent (m)"
            arg_type = FT
            default = FT(1800)
        "--ymax"
            help = "Specify maximum y extent (m)"
            arg_type = FT
            default = FT(1800)
        "--zmax"
            help = "Specify maximum z extent (m)"
            arg_type = FT
            default = FT(4000)
        "--tmax"
            help = "Specify maximum time of the simulation (s)"
            arg_type = FT
            default = FT(3600*6)
        "--moisture_model"
            help = "Moisture model to use - eq or non-eq"
            arg_type = String
            default = "nonequilibrium"
        "--precipitation_model"
            help = "precip model"
            arg_type = String
            default = "no_precipitation"

        "--tau_cond_evap"
            help = "Liquid condensation/evaporation relaxation timescale"
            arg_type = typeof(τ_cond_evap(param_set.microphys.liq)) # I think it has to be this way if you want an actual number out, rather then a method? param_set is an actual instance, but param_set has no method at the top level
            default = τ_cond_evap(param_set.microphys.liq) #FT(10) # see https://github.com/CliMA/CLIMAParameters.jl/blob/master/src/Atmos/atmos_parameters.jl
        "--tau_sub_dep"
            help = "Ice sublimation/deposition relaxation timescale"
            arg_type = typeof(τ_sub_dep(param_set.microphys.ice))
            default = τ_sub_dep(param_set.microphys.ice) #FT(10) # see https://github.com/CliMA/CLIMAParameters.jl/blob/master/src/Atmos/atmos_parameters.jl

        "--tau_cond_evap_scale"
            help = "tau_cond_evap scaling factor away from default or prescibed value"
            arg_type = FT
            default = FT(1)
        "--tau_sub_dep_scale"
            help = "τ_sub_dep scaling factor away from default or prescibed value"
            arg_type = FT
            default = FT(1)
        
        "--timestep"
            help = "Timestep for the model solver"
            arg_type = FT # let it resolve to whatever you put in (same type as others), but will default to nothing

        "--solver_type"
            help = "The solver type for the model"
            arg_type = String
            default = "imex_solver"

    end
    return ClimateMachine.init(parse_clargs = parse_clargs, custom_clargs = cfsite_args,init_driver=init_driver)

end

cl_args = parse_commandline()

# These const declarations are essential
const _τ_cond_evap       = cl_args["tau_cond_evap"]
const _τ_sub_dep         = cl_args["tau_sub_dep"]
const _τ_cond_evap_scale = cl_args["tau_cond_evap_scale"]
const _τ_sub_dep_scale   = cl_args["tau_sub_dep_scale"]

# -- I think here we're assigning a fcn method to EarthParameter Set, which is now top level in param_set, and hence it doesn't look deeper to param_set.microphys... ?
# ... not sure if this is the way but who knows... Maybe i should also overwrite the fcn in param_set.microphys... but idk if that makes a diff...
CLIMAParameters.Atmos.Microphysics.τ_cond_evap(::EarthParameterSet) = _τ_cond_evap_scale * _τ_cond_evap # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/
CLIMAParameters.Atmos.Microphysics.τ_sub_dep(::EarthParameterSet)   = _τ_sub_dep_scale   * _τ_sub_dep   # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/

# CLIMAParameters.Atmos.Microphysics.τ_cond_evap(::AbstractLiquidParameterSet) = _τ_cond_evap_scale * _τ_cond_evap # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/
# CLIMAParameters.Atmos.Microphysics.τ_sub_dep(::AbstractIceParameterSet)      = _τ_sub_dep_scale   * _τ_sub_dep   # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/

CLIMAParameters.Atmos.Microphysics.τ_cond_evap(::LiquidParameterSet) = _τ_cond_evap_scale * _τ_cond_evap # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/
CLIMAParameters.Atmos.Microphysics.τ_sub_dep(::IceParameterSet)      = _τ_sub_dep_scale   * _τ_sub_dep   # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/

# print("Running model with τ_cond_evap: ", τ_cond_evap(param_set), "\n") #  CLIMAParameters.Atmos.Microphysics.τ_cond_evap(param_set) also works
# print("Running model with τ_sub_dep: "  , τ_sub_dep(param_set)  , "\n") #  CLIMAParameters.Atmos.Microphysics.τ_sub_dep(param_set)   also works

print("Running model with τ_cond_evap: ", τ_cond_evap(param_set.microphys.liq), "\n") #  CLIMAParameters.Atmos.Microphysics.τ_cond_evap(param_set) also works
print("Running model with τ_sub_dep: "  , τ_sub_dep(param_set.microphys.ice)  , "\n") #  CLIMAParameters.Atmos.Microphysics.τ_sub_dep(param_set)   also works

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Model configuration at bottom in main()

# Citation for problem setup
"""
CMIP6 Test Dataset - cfsites
@Article{gmd-10-359-2017,
AUTHOR = {Webb, M. J. and Andrews, T. and Bodas-Salcedo, A. and Bony, S. and Bretherton, C. S. and Chadwick, R. and Chepfer, H. and Douville, H. and Good, P. and Kay, J. E. and Klein, S. A. and Marchand, R. and Medeiros, B. and Siebesma, A. P. and Skinner, C. B. and Stevens, B. and Tselioudis, G. and Tsushima, Y. and Watanabe, M.},
TITLE = {The Cloud Feedback Model Intercomparison Project (CFMIP) contribution to CMIP6},
JOURNAL = {Geoscientific Model Development},
VOLUME = {10},
YEAR = {2017},
NUMBER = {1},
PAGES = {359--384},
URL = {https://www.geosci-model-dev.net/10/359/2017/},
DOI = {10.5194/gmd-10-359-2017}
}
"""

struct GCMRelaxation{FT} <: AbstractSource
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

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Temperature tendency term, ∂T∂t
"""
    EnergyTendency <: Source

Temperature tendency for the LES configuration based on quantities 
from a GCM. Quantities are included in standard CMIP naming format. 
Tendencies included here are 

    tntha = temperature tendency due to horizontal advection
    tntva = temperature tendency due to vertical advection
    tntr = temperature tendency due to radiation fluxes
    ∂T∂z = temperature vertical gradient from GCM values
"""
struct EnergyTendency <: AbstractSource end
function atmos_source!(
    s::EnergyTendency,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    # Establish problem float-type
    FT = eltype(state)
    _grav = grav(atmos.param_set)
    # Establish vertical orientation
    k̂ = vertical_unit_vector(atmos, aux)
    _e_int_v0 = e_int_v0(atmos.param_set) # """ Specific internal energy of vapor at ``T_0`` (J/kg) """, see https://github.com/CliMA/CLIMAParameters.jl/blob/f779db6fa86f0bfce9de4ee31757fed1fe8b75ed/src/Planet/Planet.jl
    # Unpack vertical gradients
    # ∂qt∂z = dot(diffusive.gcminfo.∇hus, k̂)
    # ∂T∂z  = dot(diffusive.gcminfo.∇ta, k̂)
    # w_s   = -aux.gcminfo.wap / aux.gcminfo.ρ / _grav
    # ---------------------------------------------------------------------- new driver version | HadGEMVertical()
    ∂qt∂z = diffusive.gcminfo.∇ᵥhus
    ∂T∂z  = diffusive.gcminfo.∇ᵥta
    w_s   = aux.gcminfo.w_s
    # -----------------------------------------------------------------------------------------
    # Establish thermodynamic state
    TS = recover_thermo_state(atmos, state, aux)
    cvm = cv_m(TS)
    # Compute tendency terms
    # Temperature contribution
    # T_tendency =
        # aux.gcminfo.tntha + aux.gcminfo.tntva + aux.gcminfo.tntr + ∂T∂z * w_s
    # ---------------------------------------------------------------------- new driver version | HadGEMVertical()
    T_tendency =
        aux.gcminfo.Σtemp_tendency + ∂T∂z * w_s
    # -----------------------------------------------------------------------------------------
    # Moisture contribution
    # q_tot_tendency = aux.gcminfo.tnhusha + aux.gcminfo.tnhusva + ∂qt∂z * w_s
    # ---------------------------------------------------------------------- new driver version | HadGEMVertical()
    q_tot_tendency = aux.gcminfo.Σqt_tendency + ∂qt∂z * w_s
    # -----------------------------------------------------------------------------------------
    source.ρe += cvm * state.ρ * T_tendency
    source.ρe += _e_int_v0 * state.ρ * q_tot_tendency
    # GPU-friendly return nothing
    return nothing
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

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
struct MoistureTendency <: AbstractSource end
function atmos_source!(
    s::MoistureTendency,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    # Establish problem float-type
    FT = eltype(state)
    _grav = grav(atmos.param_set)
    k̂ = vertical_unit_vector(atmos, aux)
    # Establish vertical orientation
    # ∂qt∂z = dot(diffusive.gcminfo.∇hus, k̂)
    # w_s = -aux.gcminfo.wap / aux.gcminfo.ρ / _grav
    # ---------------------------------------------------------------------- new driver version | HadGEMVertical()
    ∂qt∂z = diffusive.gcminfo.∇ᵥhus
    w_s = aux.gcminfo.w_s # gcm vertical motion isn't great for moisture here sincei it's convective and our moisture profile will change...
    #
    # -----------------------------------------------------------------------------------------
    # Establish thermodynamic state
    TS = recover_thermo_state(atmos, state, aux)
    cvm = cv_m(TS) # gas constant
    # Compute tendency terms
    # q_tot_tendency = aux.gcminfo.tnhusha + aux.gcminfo.tnhusva + ∂qt∂z * w_s
    # ---------------------------------------------------------------------- new driver version | HadGEMVertical()
    q_tot_tendency = aux.gcminfo.Σqt_tendency + ∂qt∂z * w_s
    # -----------------------------------------------------------------------------------------
    source.moisture.ρq_tot += state.ρ * q_tot_tendency
    source.ρ               += state.ρ * q_tot_tendency
    # GPU-friendly return nothing
    return nothing
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Large-scale subsidence forcing
"""
    Subsidence <: Source 

Subsidence tendency, given a vertical velocity at the large scale, 
obtained from the GCM data. 

    wap = GCM vertical velocity [Pa s⁻¹]. Note the conversion required
"""
struct SubsidenceTendency <: AbstractSource end
function atmos_source!(
    s::SubsidenceTendency,
    atmos::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    _grav = grav(atmos.param_set)
    # Establish vertical orientation
    k̂ = vertical_unit_vector(atmos, aux)
    # Establish subsidence velocity
    # w_s = -aux.gcminfo.wap / aux.gcminfo.ρ / _grav
    # ---------------------------------------------------------------------- new driver version | HadGEMVertical()
    w_s = aux.gcminfo.w_s # gcm vertical motion isn't great for subsidence here sincei it's convective
    # -----------------------------------------------------------------------------------------

    # Compute tendency terms
    source.ρe              -= state.ρ * w_s * dot(k̂, diffusive.∇h_tot)
    source.moisture.ρq_tot -= state.ρ * w_s * dot(k̂, diffusive.moisture.∇q_tot)
    
    # GPU-friendly return nothing
    return nothing
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

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
struct LinearSponge{FT} <: AbstractSource
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
    direction,
)
    #Unpack sponge parameters
    FT = eltype(state)
    z_max = s.z_max
    z_sponge = s.z_sponge
    α_max = s.α_max
    γ = s.γ
    # Establish sponge relaxation velocity
    u_geo = SVector(aux.gcminfo.ua, aux.gcminfo.va, aux.gcminfo.w_s) # the velocity vector
    z = altitude(atmos, aux)
    # Accumulate sources
    if z_sponge <= z
        r = (z - z_sponge) / (z_max - z_sponge)
        #ZS: different sponge formulation?
        β_sponge =
            α_max .* sinpi(r / FT(2)) .^ γ

        mom_relax  = β_sponge * (state.ρu .- aux.gcminfo.ρ * u_geo)
        source.ρu -= mom_relax
        # print(source.pu, u_geo, "\n")

        # e_int_relax = .5 * (norm(state.ρu) - norm(state.ρu - mom_relax)) # KE before - KE after # dont think this works because source units are not the same...

        # also relax temperature to fixed -- 2 part maybe? how to stabilize
        # -- if we dont want to calculate the internal energy from the gcm  because it might be complicated:
        # -- -- let's just relax ta, and then calculate the temp tendency of that relaxation and the associated energy contribution
        # -- -- moisture then gets adjusted to remain in equilibrium but it is probably undersaturated and is fine... (I hope...)
        
        # implementation (hard bc there's no temp variable it's calculated from the thermodynamic state -- idk seems janky)
        ts = recover_thermo_state(atmos, state, aux)       
        cvm = cv_m(ts)                                  
        T  = air_temperature(ts) 
        # calculate temperature tendency from this, and then calculate forcing on ρe, but no way to guarantee if affects temp and not moisture?
        T_tendency  =  (T              .-  aux.gcminfo.ta ) 
        source.ρe  -= β_sponge * (cvm * state.ρ * T_tendency) #.- e_int_relax # gauranteed to work I suppose exluding moisture...


        # alternatively, if we could could the internal energy from the reference (gcm) state:
        # -- maybe relax internal energy ρe to its gcm value, and then relax total moisture source.moisture.ρq_tot or ta to its value such that the thermodynamic state is fixed
        # -- -- fixing both might have the unintended effect of being out of equilibrium so we'll try out best to remain in eq this way by doing one or the other
        source.moisture.ρq_tot -= β_sponge * (state.moisture.ρq_tot .-  aux.gcminfo.ρq_tot) # relaxes moisture maybe?      see , moisture tendency for energy

        # I think what happens is T can drift the wrong way but the correction to internal energy somehow dosnn't make it to energy and adjust the state? idk it's weird
        # print(z, " | ", state.ρe ,  " | ", aux.gcminfo.ρe, " | ", -mom_relax, " | ", "\n")
        # print(-β_sponge * (cvm * state.ρ * T_tendency), " | ", source.ρe ," | ", cvm * state.ρ * aux.gcminfo.Σtemp_tendency, " | ",  -β_sponge * (state.ρe   .-  aux.gcminfo.ρe) ,    "\n" )  
        # print("------", "\n")

        # source.ρe              -= β_sponge * (state.ρe              .-  aux.gcminfo.ρe    ) #.-  e_int_relax # relaxe total energy #maybe this is redundant if you do moisture and velocity
        
        # Not sure how to account for moisture changes but yah is ok I suppose... idk if is better or wors ethan 
    
        # I'm not sure if this helps but it might not hurt with fixing ta in the model towards what it is in the gcm...
        # source.ρ               -= β_sponge * (state.ρ               .-  aux.gcminfo.ρ ) # relax total density ... might break the model (actually if you don't do this seems T will still drift, idk)
        # maybe density relaxation creates our gravity waves> try without it and see if can stop temp from drifting anyway? idk

        # GPU-friendly return nothing

    end


    # GPU-friendly return nothing
    return nothing
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# We first specify the NetCDF file from which we wish to read our 
# GCM values.
# Utility function to read and store variables directly from the 
# NetCDF file
"""
    str2var(str::String, var::Any)

Helper function allowing variables read in from the GCM file
to be made available to the LES simulation.
"""
function str2var(str::String, var::Any)
    str = Symbol(str)
    @eval(($str) = ($var))
end

# Define the get_gcm_info function
"""
    get_gcm_info(...)

For a specific global site, establish and store the GCM state
for each available vertical level. refers to the integer
index of the specific global site that we are interested in.
"""
function get_gcm_info(; # should be overwritten by default using command_line arguments
                      data_path = "",  
                      model     = "",
                      exper     = "",
                      rip       = "",
                      sites     = "all", 
                      years     = "all",
                      months    = "all",
                      days      = "all",
                      hours     = "all",
                      minutes   = "all",
                      seconds   = "all"
                      )

    @printf("--------------------------------------------------\n")
    @info @sprintf("""\n
     Experiment: GCM-driven LES(ClimateMachine)
     """)

    @printf("\n")
    print("LES = " * model * "-"* exper * "-" *rip * "  |  site: " * string(sites) * " years: " * string(years) * " months: " * string(months) * " days: " * string(days) * " hours: " * string(hours) * " minutes: " * string(minutes) * " seconds: " * string(seconds) * "\n")
    @printf("--------------------------------------------------\n")
    filepath = data_path * "/" * model * "/" * exper * "/"
        # "/central/groups/esm/zhaoyi/GCMForcedLES/forcing/clima/" *
        # forcingfile *
        # ".nc"
    print(filepath * "\n")
    # filenames = filter(contains(r".nc"), readdir(filepath,join=true)) # works in 1.5 but deprecated before?
    filenames = filter(x->occursin(".nc",x), readdir(filepath,join=true)) # julia 1.4.2
    req_varnames = (
        "zg",
        "ta",
        "hus",
        "ua",
        "va",
        "pfull",
        "tntha",
        "tntva",
        "tntr",
        "tnhusha",
        "tnhusva",
        "wap",
        "hfls",
        "hfss",
	    "ts",
        "alpha",
        "cli", # taken to help initialize the model... (noneq seems to not like initializing w/o it wit imex_solver)
        "clw"
    )
    # Load NETCDF dataset (HadGEM2-A information)
    # Load the NCDataset (currently we assume all time-stamps are 
    # in the same NCData file). We store this information in `data`. 

    # NCDataset below chokes on merging empty w/ nonempty files so we try to fix that
    # data     = Array{Any}(undef,length(filenames))
    times    = Array{Any}(undef,length(filenames)) # loads times incorrectly (seems to copy from first file)
    indices  = collect(1:length(filenames))

    for (index, file) in enumerate(filenames)
        data = NCDataset(file)
        if length(data["time"]) == 0
            indices = filter!(e->e≠index,indices) # drop those empty files
        else
            times[index] = data["time"] # save the right times
        end
    end
    
    times = times[indices]
    filenames = filenames[indices]
    time = cat(times...;dims=1) # unpacks using ellipsis

    data = NCDataset(filenames, aggdim = "time", deferopen = false)  # loads times incorrectly (seems to copy from first file)

    if years === "all"
        years = Dates.year.(time) # get all the years
        years = Array(minimum(years):1:maximum(years)) # get an individula vector
    end  
    #
    if months === "all" # double or triple equals? also is there an `is` operator?
        # months =  Array(1:1:12)
        months = Dates.month.(time) # get all the years
        months = Array(minimum(month):1:maximum(month)) # get an individula vector
    end
    #
    if days === "all" # double or triple equals? also is there an `is` operator?
        # dates =  Array(1:1:31)
        days = Dates.day.(time) # get all the years
        days = Array(minimum(days):1:maximum(days)) # get an individula vector
    end
    #
    if hours === "all" # double or triple equals? also is there an `is` operator?
        # hours =  Array(0:1:23)
        hours = Dates.hour.(time) # get all the years
        hours = Array(minimum(hours):1:maximum(hours)) # get an individula vector
    end
    #
    if minutes === "all" # double or triple equals? also is there an `is` operator?
        # minutes =  Array(0:1:59)
        minutes = Dates.minute.(time) # get all the years
        minutes = Array(minimum(minutes):1:maximum(minutes)) # get an individula vector
    end
    #
    if seconds === "all" # double or triple equals? also is there an `is` operator?
        # minutes =  Array(0:1:59)
        seconds = Dates.second.(time) # get all the years
        seconds = Array(minimum(seconds):1:maximum(seconds)) # get an individula vector
    end



    time_mask =  (Dates.year.(time)   .∈  Ref(years))   .& 
                 (Dates.month.(time)  .∈  Ref(months))  .& 
                 (Dates.day.(time)    .∈  Ref(days))    .& 
                 (Dates.hour.(time)   .∈  Ref(hours))   .& 
                 (Dates.minute.(time) .∈  Ref(minutes)) .& 
                 (Dates.second.(time) .∈  Ref(seconds)) 


    if sites === "all"
        sites = data["site"][:]
    end
    site_mask = data["site"] .∈ Ref(sites)

    print("time_mask sum: " , string(sum(time_mask)),   " | site_mask_sum: ", string(sum(time_mask)),"\n")

    # print(data)

    # To assist the user / inform them of the data processing step
    # we print out some useful information, such as groupnames 
    # and a list of available variables
    @printf("Storing information for group %s ...\n", sites)
    for (varname, var) in data #.group[groupid]
        for reqvar in req_varnames
            if reqvar == varname
                print("handling " * varname,"\n")
                if varname == "hfls" || varname == "hfss" || varname == "ts" # surface properties, have no lev as dim 1
                    var = var[:,:]
                    var = var[site_mask,time_mask]
                    var = Statistics.mean(var, dims = [1,2])[1] # should work
                else
                    var = var[:,:,:] # Loading in advance makes it much faster
                    # print(size(var),"\n")
                    var = var[:,site_mask,time_mask] # seems slow for some reason, also check order
                    # print(size(var),"\n")
                    # Get average over time dimension
                    var = Statistics.mean(var, dims = [2, 3])
                end
                # print(size(var),"\n")
                # Assign the variable values to the appropriate converted string
                str2var(varname, var)
            end
        end
        # Store key variables
    end
    @printf("Complete\n")
    @printf("--------------------------------------------------\n")
    # @printf("Group data storage complete\n")

    
    return (
        zg,
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
        cli,
        clw
    )

end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Initialise the CFSite experiment :D! (where is this called with these args?)
const seed = MersenneTwister(0)
function init_cfsites!(problem, bl, state, aux, localgeo, t, spl)
    FT = eltype(state)
    (x, y, z) = localgeo.coord # added from new branch
    _grav = grav(bl.param_set)

    # Unpack splines, interpolate to z coordinate at 
    # present grid index. (Functions are all pointwise)

    # print("spl_ta: "  ,spl.spl_ta , " ", z, "\n")
    # print("spl_hus: " ,spl.spl_hus, " ", z, "\n")


    ta      = FT(spl.spl_ta(z))
    hus     = FT(spl.spl_hus(z))
    ua      = FT(spl.spl_ua(z))
    va      = FT(spl.spl_va(z))
    pfull   = FT(spl.spl_pfull(z))
    tntha   = FT(spl.spl_tntha(z))
    tntva   = FT(spl.spl_tntva(z))
    tntr    = FT(spl.spl_tntr(z))
    tnhusha = FT(spl.spl_tnhusha(z))
    tnhusva = FT(spl.spl_tnhusva(z))
    wap     = FT(spl.spl_wap(z))
    ρ_gcm   = FT(1 / spl.spl_alpha(z))

    cli     = FT(spl.spl_cli(z))
    clw     = FT(spl.spl_clw(z))

    w_s   = -wap / ρ_gcm / _grav 

    # print("ta: " ,ta, " ","hus: " ,hus," ", z, "\n")
    include_vertical_velocity = FT(false)

    # Compute field properties based on interpolated data
    ρ     =     air_density(bl.param_set, ta, pfull, PhasePartition(hus))
    e_int = internal_energy(bl.param_set, ta       , PhasePartition(hus))
    e_kin = (ua^2 + va^2 + (w_s*include_vertical_velocity)^2) / 2
    e_pot = _grav * z

    uv    = SVector(ua, va, w_s*include_vertical_velocity) # initialize the velocity vector
    e_tot = e_kin + e_pot + e_int
    q_tot = hus + clw + cli
    # Assignment of state variables
    state.ρ  = ρ
    state.ρu = ρ * uv # we added in  vertical velocity but could also not init wit it
    state.ρe = ρ * e_tot
    # state.moisture.ρq_tot = ρ * hus
    state.moisture.ρq_tot = ρ * q_tot # initialize wit all the moisture bruh


    if bl.moisture isa NonEquilMoist # added from bomex, intialize to 0 if nonequilibrium (can try initializing to model also one day)
        # state.moisture.ρq_liq = FT(0)
        # state.moisture.ρq_ice = FT(0)
        state.moisture.ρq_liq = clw # tryn this out
        state.moisture.ρq_ice = cli

    if bl.precipitation isa Rain
        state.precipitation.ρq_rai = FT(0)
    end 



        # q_tot = hus
        # q_pt = PhasePartition(q_tot, cli, clw)      # (q_tot::Real[, q_liq::Real[, q_ice::Real]]), not sure if just hus or ρq_tot

    
        # sol = ClimateMachine.Thermodynamics.PhaseNonEquil_ρTq( # in Thermodynamics/states, see https://clima.github.io/ClimateMachine.jl/latest/APIs/Common/Thermodynamics/#ClimateMachine.Thermodynamics.PhaseNonEquil_%CF%81Tq)
            # param_set, ρ, ta, q_pt)  # seems to adjust the state w/o doing anything, idk...

        # returns PhaseNonEquil{FT, typeof(param_set)}(param_set, e_int, ρ, q_pt)


        # sol = saturation_adjustment( param_set, e_int, 
        #     param_set    , # param_set::APS,
        #     e_int        , # e_int::FT,
        #     ρ            , # ρ::FT,
        #     q_tot        , # q_tot::FT,
        #     PhaseNonEquil, # phase_type::Type{<:PhaseEquil},    < not sure what to choose here >
        #     5            , # maxiter::Int, < maxiter = 5 >
        #     FT(0.1)      , # temperature_tol::FT, < tolerance = FT(0.1) >
        # ) 

        # print(sol)


    end
    # print("state before: ", state.moisture.ρq_tot," ",ρ, " ",  hus," ", state.moisture, "\n")
    # print("-----------","\n")
    if z <= FT(400)
        state.ρe += rand(seed) * FT(1 / 100) * (state.ρe)
        state.moisture.ρq_tot +=
            rand(seed) * FT(1 / 100) * (state.moisture.ρq_tot)
    end
    # print("state after: ", state.moisture.ρq_tot," ",ρ, " ",  hus, " ", state.moisture, "\n")

    # Assign and store the ref variable for sources
    # aux.gcminfo.ρ       = ρ_gcm
    # aux.gcminfo.pfull   = pfull
    # aux.gcminfo.ta      = ta
    # aux.gcminfo.hus     = hus
    # aux.gcminfo.tntha   = tntha
    # aux.gcminfo.tntva   = tntva
    # aux.gcminfo.ua      = ua
    # aux.gcminfo.va      = va
    # aux.gcminfo.tntr    = tntr
    # aux.gcminfo.tnhusha = tnhusha
    # aux.gcminfo.tnhusva = tnhusva
    # aux.gcminfo.wap     = wap

    # ---------------------------------------------------------------------- new driver version | HadGEMVertical()
    # Assign and store the ref variable for sources

    tntr_value = FT(-1/86400) # K/day * 1d/86400s


    # since we defined this in terms of the splines, we should use the ρ not the ρ_gcm later
    aux.gcminfo.ta             = ta # 
    aux.gcminfo.hus            = hus
    aux.gcminfo.Σtemp_tendency = (tntha + tntva + tntr) * 0 + tntr * tntr_value/tntr # ignore the temperature fluxes from advection since they might be large and random, radiation pretty irrelevant
    aux.gcminfo.ua             = ua
    aux.gcminfo.va             = va
    aux.gcminfo.Σqt_tendency   = (tnhusha + tnhusva) * 0 # Added to try to reign in the madness, no fluxes of moisture from environment
    aux.gcminfo.w_s            = w_s * include_vertical_velocity  # this is bad for forcing since it can be large in covective states, kill

    # ρ     =     air_density(bl.param_set, aux.gcm.ta, pfull, PhasePartition(hus))
    e_int_gcm = internal_energy(bl.param_set, aux.gcminfo.ta       , PhasePartition(aux.gcminfo.hus))
    e_pot_gcm = _grav * z
    e_kin_gcm = (aux.gcminfo.ua^2 + aux.gcminfo.va^2 + aux.gcminfo.w_s^2) / 2
    e_tot_gcm = e_int_gcm + e_pot_gcm + e_kin_gcm
    # --> added to help with sponge boundary conditions and keeping track of this stuff in the gcm state... not sure if e_int will work tho...
    # -- --> not sure if to use ρ_gcm or ρ for these, i guess it depends on what partition u want... I'll go with the un adjusted ones since might vibe wit ta better?
    # print("e_int is: ", e_int, ", FT is: ", string(FT), ", Resolution is: ", FT(e_int + ua*0) , ", FT(ua) is: ", FT(ua), "\n")
    # aux.gcminfo.e_int           = FT(e_int + ua*0) # test if it's a type thing and addding ua*a can help
    # aux.gcminfo.e_kin           = FT(e_kin)
    # aux.gcminfo.e_pot           = FT(e_pot)
    
    # use ρ not ρ_gcm so the relaxation works as intended towards the initial state w/o instabilities....
    aux.gcminfo.ρe                = ρ * e_tot_gcm
    aux.gcminfo.ρ                 = ρ
    aux.gcminfo.ρq_tot            = ρ * q_tot
    # ------------------------------------------------------------------------------------------

    return nothing
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

function config_cfsites( # returns config 
    FT,
    N,
    resolution,
    xmax,
    ymax,
    zmax,
    hfls,
    hfss,
    ts,
    sites,
    moisture_model,
    precipitation_model,
    solver_type,
)


    test_scale = FT(2.0)

    # Boundary Conditions
    u_star = FT(0.28) # Friction velocity (explanation from bomex_model)

    print("hfls: ",hfls,"\n")
    print("hfss: ",hfss,"\n")

    problem = AtmosProblem(
        boundaryconditions = (
            AtmosBC(
                momentum = Impenetrable(DragLaw(
                    (state, aux, t, normPu_int) -> (u_star / normPu_int)^2,
                )),
                energy = PrescribedEnergyFlux((state, aux, t) -> (abs(hfls) + abs(hfss)) * test_scale),
                moisture = PrescribedMoistureFlux(
                    (state, aux, t) ->
                        abs(hfls)*test_scale / latent_heat_vapor(param_set, ts),
                ),
            ),
            AtmosBC(),
        ),
        init_state_prognostic = init_cfsites!,
    )


    source_default = (
        Gravity(),
        LinearSponge{FT}(zmax, zmax * 0.85, 1, 4),
        MoistureTendency(),
        EnergyTendency(),
        SubsidenceTendency(),
    )
    
    # moisture = EquilMoist{FT}(; maxiter = 5, tolerance = FT(2)),
    # FROM BOMEX
    if moisture_model == "equilibrium"
        # print("Using equilibrium moisture model..." * "\n")
        source = source_default
        moisture = EquilMoist{FT}(; maxiter = 5, tolerance = FT(0.1))
    elseif moisture_model == "nonequilibrium"
        # print("Using nonequilibrium moisture model..." * "\n")
        # source = (source_default..., CreateClouds())
        # source = (source_default..., CreateClouds(),RemovePrecipitation(true)) # testing w/ precip?
        source = (source_default..., CreateClouds()) # testing w/ precip?

        moisture = NonEquilMoist()
    else
        @warn @sprintf(""" %s: unrecognized moisture_model in source terms, using the defaults""",
            moisture_model,
        )
        moisture = EquilMoist{FT}(; maxiter = 5, tolerance = FT(0.1))
        source = source_default
    end

    print("Using " , moisture_model , " moisture model: " , moisture , "\n")


    # precipitation model and its sources
    if precipitation_model == "no_precipitation"
        precipitation = NoPrecipitation()
    elseif precipitation_model == "rain"
        source = (source..., Rain_1M())
        precipitation = Rain()
    elseif precipitation_model == "remove_precipitation"
        source = (source..., RemovePrecipitation(true))
        precipitation = NoPrecipitation()
    else
        @warn @sprintf(
            """%s: unrecognized precipitation_model in source terms, using the defaults""",
            precipitation_model,
        )
        precipitation = NoPrecipitation()
    end

    print("Using " , precipitation_model , " moisture model: " , precipitation , "\n")


    model = AtmosModel{FT}(
        AtmosLESConfigType,
        param_set;
        problem = problem,
        turbulence = Vreman{FT}(0.23),
        source = source,
        moisture = moisture,
        precipitation = precipitation,
        #ZS: hyperdiffusion?
        #hyperdiffusion = DryBiharmonic{FT}(12*3600),
        # gcminfo = HadGEM(), #????
        # gcminfo = HadGEMVertical(), # in new branch is HadGEMVertical() but where is this defined\
        gcminfo = CMIP_cfsite_Vertical(), # my own adapted version to keep the variables I might need
    )

    # # Timestepper options
    # # ex_solver   = ClimateMachine.ExplicitSolverType() 
    # # imex_solver = ClimateMachine.IMEXSolverType()
    # if     solver_type == "ex_solver"
    #     solver_type = ClimateMachine.ExplicitSolverType() # Explicit Solver
    # elseif solver_type == "imex_solver"
    #     solver_type = ClimateMachine.IMEXSolverType()     # IMEX Solver Type
    # elseif solver_type == "mrrk_solver"
    #     solver_type = ClimateMachine.MultirateSolverType( # Multirate Explicit Solver
    #     fast_model = AtmosAcousticGravityLinearModel,
    #     slow_method = LSRK144NiegemannDiehlBusch,
    #     fast_method = LSRK144NiegemannDiehlBusch,
    #     timestep_ratio = 10)
    # elseif solver_type == "LSRK144NiegemannDiehlBusch"
    #     solver_type = ClimateMachine.ExplicitSolverType( # why no semicolon for keyword arg here?
    #     solver_method = LSRK144NiegemannDiehlBusch,
    #     )
    # else
    #     @warn @sprintf(""" %s: unrecognized solver type, using the defaults""",
    #     solver_type = ClimateMachine.IMEXSolverType() 
    # )
    # end
    # print("Using ", solver_type, " solver type...","\n" )

    # Configuration
    config = ClimateMachine.AtmosLESConfiguration(
        "sites_"*string(sites), # is just a name, see below (is used in file output so sites_# works)
        # forcingfile * "_$groupid", # is just a name https://github.com/CliMA/ClimateMachine.jl/blob/553d0568f21e376fa0ec35d83a5edacd327ebf41/src/Driver/driver_configs.jl
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        param_set,
        init_cfsites!,
        #ZS: multi-rate? -- maybe the mrrk goes here?
        solver_type = solver_type,
        model = model,
    )
    return config
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

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

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

function main()

    # Working precision
    FT = Float64

    cl_args = parse_commandline(;parse_clargs = true, init_driver=true)

    # cl_args =
    #     ClimateMachine.init(parse_clargs = true, custom_clargs = cfsite_args)
    # groupid = cl_args["group_id"]

    data_path      = cl_args["data_path"]
    model          = cl_args["model"]
    exper          = cl_args["exper"]
    rip            = cl_args["rip"]
    years          = cl_args["years"]
    months         = cl_args["months"]
    days           = cl_args["days"]
    hours          = cl_args["hours"]
    minutes        = cl_args["minutes"]
    seconds        = cl_args["seconds"]
    sites          = cl_args["sites"]

    if years === "all"
        nothing
    else
        years = eval(Meta.parse(years))
    end

    if months === "all"
        nothing
    else
        months = eval(Meta.parse(months))
    end

    if days === "all"
        nothing
    else
        days = eval(Meta.parse(days))
    end

    if hours === "all"
        nothing
    else
        hours = eval(Meta.parse(hours))
    end

    if minutes === "all"
        nothing
    else
        minutes = eval(Meta.parse(minutes))
    end

    if seconds === "all"
        nothing
    else
        seconds = eval(Meta.parse(seconds))
    end

    if sites === "all"
        nothing
    else
        sites = eval(Meta.parse(sites))
    end


    # Domain resolutions
    Δh         = cl_args["delta_h"]
    Δv         = cl_args["delta_v"]
    resolution = (Δh, Δh, Δv)
    # Domain extents
    xmax       = cl_args["xmax"]
    ymax       = cl_args["ymax"]
    zmax       = cl_args["zmax"]
    # Simulation time
    t0         = FT(0)
    timeend    = cl_args["tmax"]
    timestep   = cl_args["timestep"]


    # Microphysics :: Moisture, and relaxation parameters
    moisture_model      = cl_args["moisture_model"]
    precipitation_model = cl_args["precipitation_model"] #"rain"

    solver_type         = cl_args["solver_type"]

    print("solver_type is: " , solver_type, "\n")

    _τ_cond_evap        = cl_args["tau_cond_evap"]
    _τ_sub_dep          = cl_args["tau_sub_dep"]
    _τ_cond_evap_scale  = cl_args["tau_cond_evap_scale"]
    _τ_sub_dep_scale    = cl_args["tau_sub_dep_scale"]


    # ---- Moved to top since complains about making a global definition not in global scope... Has to be written this way I think tho bc is a multiple dispatch overload of this fcn for this type 
    # CLIMAParameters.Atmos.Microphysics.τ_cond_evap(::EarthParameterSet) = _τ_cond_evap_scale * _τ_cond_evap # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/
    # CLIMAParameters.Atmos.Microphysics.τ_sub_dep(::EarthParameterSet)   = _τ_sub_dep_scale   * _τ_sub_dep   # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/


    # CLIMAParameters.Atmos.Microphysics.τ_cond_evap(::AbstractLiquidParameterSet) = _τ_cond_evap_scale * _τ_cond_evap # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/
    # CLIMAParameters.Atmos.Microphysics.τ_sub_dep(::AbstractIceParameterSet)      = _τ_sub_dep_scale   * _τ_sub_dep   # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/
    
    # CLIMAParameters.Atmos.Microphysics.τ_cond_evap(::LiquidParameterSet) = _τ_cond_evap_scale * _τ_cond_evap # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/
    # CLIMAParameters.Atmos.Microphysics.τ_sub_dep(::IceParameterSet)      = _τ_sub_dep_scale   * _τ_sub_dep   # for some reason something like τ_cond_evap(param_set) fails (recursive memory error?), shrug, anyway see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/

    print("Running model with τ_cond_evap: ", τ_cond_evap(param_set), "\n") #  CLIMAParameters.Atmos.Microphysics.τ_cond_evap(param_set) also works
    print("Running model with τ_sub_dep: "  , τ_sub_dep(param_set)  , "\n") #  CLIMAParameters.Atmos.Microphysics.τ_sub_dep(param_set)   also works

    print("Running model with τ_cond_evap: ", τ_cond_evap(param_set.microphys.liq), "\n") #  CLIMAParameters.Atmos.Microphysics.τ_cond_evap(param_set) also works
    print("Running model with τ_sub_dep: "  , τ_sub_dep(param_set.microphys.ice)  , "\n") #  CLIMAParameters.Atmos.Microphysics.τ_sub_dep(param_set)   also works



    # or if default was nothing, perhaps:
    # if !isnothing(tau_cond_evap):
    #         CLIMAParameters.Atmos.Microphysics.τ_cond_evap(::EarthParameterSet) = τ_cond_evap_scale * τ_cond_evap(param_set.microphys.liq) # see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/
    # if !isnothing(tau_cond_evap):
    #         CLIMAParameters.Atmos.Microphysics.τ_sub_dep(::EarthParameterSet)   = τ_cond_evap_scale * τ_sub_dep(param_set.microphys.ice) # see https://juliahub.com/docs/CLIMAParameters/B1Qj2/0.1.8/


    # if timestep === "timestep"
    #     nothing
    # else
    #     timestep = eval(Meta.parse(timestep))
    # end


    # DG polynomial order
    N = 4
    # Courant number
    CFL = FT(0.8)
    # CFL = FT(0.2)
    if     solver_type == "ex_solver"
        CFL = FT(0.8)
    elseif solver_type == "imex_solver"
        CFL = FT(0.2)
    elseif solver_type == "mrrk_solver"
        CFL = FT(0.8)
    elseif solver_type == "LSRK144NiegemannDiehlBusch"
        CFL = FT(0.8)
    else
        @warn @sprintf(""" %s: unrecognized solver type, using the defaults""",
        CFL = FT(0.8) 
    )
    end
    print("Using ", solver_type, " solver type...","\n" )


    # Execute the get_gcm_info function
    (   z,
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
        cli,
        clw
    ) = get_gcm_info(;
                    data_path=data_path,
                    model=model,
                    exper=exper,
                    rip=rip,
                    years=years,
                    months=months,
                    days=days,
                    hours=hours,
                    minutes=minutes,
                    seconds=seconds,
                    sites=sites
                    ) # call on init fcn


    z= z[:] # i think more consistent than dropdims
    # print("len of z is: ", length(z), " z: ", z,"\n")
    # print("ta: "  ,vec(ta), "\n")
    # print("hus: " ,vec(hus), "\n")

    ts = ts + FT(0)
    t_trans = 220
    gamma_s   = 7.5
    gamma_a   = -1
    z_trans = (ts - t_trans) / (gamma_s/1000)

    ta = (z .<= z_trans) .* (ts      .- (gamma_s *  z./1000)) + (z .> z_trans) .* (t_trans .- (gamma_a *  (z.-z_trans)./1000))
    # such a change might warrant the recalculation of density, but maybe ok since it's recalc from pfull and we discard the gcm value...
    # asuming the changes in pressure arent huge and buff themselves out with time? we don't relax rho so hopefullly not too many grav waves

    # ta = max.( ts .- (8 * z ./ 1000), 220 ) # adopt a profile somewhere between dry n moist adiabatic :P


    splines = (
        spl_ta      = Spline1D(z, vec(ta)), # i think [:], vec more consistent than view with unknown trailing singletons, not sure which is better
        spl_pfull   = Spline1D(z, vec(pfull)),
        spl_ua      = Spline1D(z, vec(ua)),
        spl_va      = Spline1D(z, vec(va)),
        spl_hus     = Spline1D(z, vec(hus)),
        spl_tntha   = Spline1D(z, vec(tntha)),
        spl_tntva   = Spline1D(z, vec(tntva)),
        spl_tntr    = Spline1D(z, vec(tntr)),
        spl_tnhusha = Spline1D(z, vec(tnhusha)),
        spl_tnhusva = Spline1D(z, vec(tnhusva)),
        spl_wap     = Spline1D(z, vec(wap)),
        spl_alpha   = Spline1D(z, vec(alpha)),
        spl_cli     = Spline1D(z, vec(cli)),
        spl_clw     = Spline1D(z, vec(clw)),
    )



    # Timestepper options
    # ex_solver   = ClimateMachine.ExplicitSolverType() 
    # imex_solver = ClimateMachine.IMEXSolverType()
    if     solver_type == "ex_solver"
        solver_type   = ClimateMachine.ExplicitSolverType() # Explicit Solver
        CFL_direction = nothing
    elseif solver_type == "imex_solver"
        solver_type   = ClimateMachine.IMEXSolverType()     # IMEX Solver Type
        CFL_direction = HorizontalDirection()
        #
        #
    elseif solver_type == "mrrk_solver"
        solver_type = ClimateMachine.MultirateSolverType( # Multirate Explicit Solver
        fast_model = AtmosAcousticGravityLinearModel,
        slow_method = LSRK144NiegemannDiehlBusch,
        fast_method = LSRK144NiegemannDiehlBusch,
        timestep_ratio = 10)
        CFL_direction = nothing # unsure
    elseif solver_type == "LSRK144NiegemannDiehlBusch"
        solver_type = ClimateMachine.ExplicitSolverType( # why no semicolon for keyword arg here?
        solver_method = LSRK144NiegemannDiehlBusch,
        )
        CFL_direction = nothing
    else
        @warn @sprintf(""" %s: unrecognized solver type, using the defaults""",
        solver_type = ClimateMachine.IMEXSolverType() 
    )
    end
    print("Using ", solver_type, " solver type...","\n" )


    # Set up driver configuration
    driver_config = config_cfsites( # is config from config_cfsites, which includes init_cfsites etc
        FT,
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        hfls,
        hfss,
        ts,
        sites,
        moisture_model,
        precipitation_model,
        solver_type
    )
    # Set up solver configuration
    # ode_dt = nothing,
    if !isnothing(CFL_direction)
        solver_config = ClimateMachine.SolverConfiguration(
            t0,
            timeend,
            driver_config,
            splines;
            ode_dt = timestep, # added in for config purposes since default seems to evaluate to nan for some reason (is from config setup, default is nothing)
            init_on_cpu = true,
            Courant_number = CFL,
            CFL_direction = CFL_direction # added from as/gcmforcing-cfsite as per zhaoyi instruction
        )
    else # use default (is there a better way to do this and convert nothing to the default?)
        solver_config = ClimateMachine.SolverConfiguration(
            t0,
            timeend,
            driver_config,
            splines;
            ode_dt = timestep, # added in for config purposes since default seems to evaluate to nan for some reason (is from config setup, default is nothing)
            init_on_cpu = true,
            Courant_number = CFL
        )
    end

    #-----------------bomex-les---------------------------------------#
    
    if moisture_model == "equilibrium"
        filter_vars = ("moisture.ρq_tot",)
    elseif moisture_model == "nonequilibrium"
        filter_vars = ("moisture.ρq_tot", "moisture.ρq_liq", "moisture.ρq_ice")
    end

    # cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
    #     Filters.apply!(
    #         solver_config.Q,
    #         filter_vars,
    #         solver_config.dg.grid,
    #         TMARFilter(),
    #     )
    #     nothing
    # end

    check_cons = (
        ClimateMachine.ConservationCheck("ρ", "3000steps", FT(0.0001)),
        ClimateMachine.ConservationCheck("ρe", "3000steps", FT(0.0025)),
    )

    #-----------------bomex-les---------------------------------------#


    # Set up diagnostic configuration
    dgn_config = config_diagnostics(driver_config)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            filter_vars, # COPIED FROM BOMEX LES (SEE ABOVE, replaces line below)
            # ("moisture.ρq_tot",),
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
