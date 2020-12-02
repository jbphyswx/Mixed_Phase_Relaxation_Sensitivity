# # Southern Ocean Mixed Phase cfSites simulations (Adapted from Dry Rayleigh Benard)

# ## Problem description
#
# 1) We wanna have LES simulations in the Southern Ocean (matching up probably w/ cfSite 92, 82, 99, maybe 78? Gotta ask Tapio )
#
# 2) Boundaries - `Sides` : Periodic? I guess...depends on I.C.s I suppose?
#                 `Top`   : Prescribed temperature, no-slip (probably fine? Assuming top is stratospheric and above convective tops)
#                 `Bottom`: Initialized temperature? (Some SST Configuration in-line w/ SST?) Should we keep this variable for more reasonable convection or fixed?
#
# 3) Domain - Probably same as the cfSites? (Tapio said start around  40km[horizontal] x 40km[horizontal] x 15km[vertical])
#
# 4) Timeend - Tapio said start ~1hr
#
# 5) Mesh Aspect Ratio (Effective Resolution) 1:1?
#
# 6) IC from cfsites? Can't seem to download data from portal... (https://esgf-node.llnl.gov/search/cmip6/)
#
# 7) No idea what to do w/ defaults? copied from 
#
# 8) Default settings can be found in src/Driver/Configurations.jl

# ## Loading code
using Distributions
using Random
using StaticArrays
using Test
using DocStringExtensions
using Printf

using ClimateMachine
ClimateMachine.init()
using ClimateMachine.Atmos #(assuming I don't need ocean? Just fix surface to some BC based on saturation or something? idk...)
using ClimateMachine.Orientations
using ClimateMachine.ConfigTypes
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.Diagnostics
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.Thermodynamics:
    TemperatureSHumEquil_given_pressure, internal_energy # something for microphsyics here????
using ClimateMachine.TurbulenceClosures
using ClimateMachine.VariableTemplates

# Include Microphysics...
using ClimateMachine.Microphysics # https://clima.github.io/ClimateMachine.jl/latest/Theory/Atmos/Microphysics/ alleges you can do this, using ClimateMachine.Atmos.Parameterizations.CloudPhysics.Microphysics #??? This is the structure that's in .../src/ClimateMachine.jl for include but idk

using CLIMAParameters
using CLIMAParameters.Planet: R_d, cp_d, cv_d, grav, MSLP
struct EarthParameterSet <: AbstractEarthParameterSet end

using CLIMAParameters.Atmos.Microphysics
struct EarthParameterSet <: AbstractEarthParameterSet end
struct LiquidParameterSet <: AbstractLiquidParameterSet end
struct RainParameterSet <: AbstractRainParameterSet end
struct IceParameterSet <: AbstractIceParameterSet end
struct SnowParameterSet <: AbstractSnowParameterSet end


const param_set = EarthParameterSet()

const liquid_param_set = LiquidParameterSet()
const rain_param_set = RainParameterSet()
const ice_param_set = IceParameterSet()
const snow_param_set = SnowParameterSet()

# Convenience struct for sharing data between kernels
struct AtmosLESConfigType{FT} # LES config type from /src/Common/ConfigTypes? Unsure if this is the right way to do this...
    xmin::FT
    ymin::FT
    zmin::FT
    xmax::FT
    ymax::FT
    zmax::FT
    T_bot::FT
    T_lapse::FT
    T_top::FT
end

# Define initial condition kernel
function init_problem!(bl, state, aux, (x, y, z), t) # somehow we gotta initialize an unstable convective profile here.....
    dc = bl.data_config
    FT = eltype(state)

    _R_d::FT = R_d(bl.param_set)
    _cp_d::FT = cp_d(bl.param_set) # add cp_m here for water?
    _grav::FT = grav(bl.param_set)
    _cv_d::FT = cv_d(bl.param_set) # add cp_m here for water?
    _MSLP::FT = MSLP(bl.param_set)

    γ::FT = _cp_d / _cv_d # we should need another for water here? 
    # δT =
    #     sinpi(6 * z / (dc.zmax - dc.zmin)) *
    #     cospi(6 * z / (dc.zmax - dc.zmin)) + rand() # is this just an almost sinusoidlal perterbation? We shouldn't need this
    # δw =
    #     sinpi(6 * z / (dc.zmax - dc.zmin)) *
    #     cospi(6 * z / (dc.zmax - dc.zmin)) + rand()
    # ΔT = _grav / _cv_d * z + δT
    
    T = dc.T_bot - ΔT # get this from some input file?
    P = _MSLP * (T / dc.T_bot)^(_grav / _R_d / dc.T_lapse)
    ρ = P / (_R_d * T)

    q_tot = FT(0)
    e_pot = gravitational_potential(bl.orientation, aux)
    ts = TemperatureSHumEquil_given_pressure(bl.param_set, T, P, q_tot)          # Should this be set with output from some file? maybe take all the mean values from an ensemble representation of a cfSite or?

    ρu, ρv, ρw = FT(0), FT(0), ρ * δw

    e_int = internal_energy(ts)
    e_kin = FT(1 / 2) * δw^2

    ρe_tot = ρ * (e_int + e_pot + e_kin)
    state.ρ = ρ
    state.ρu = SVector(ρu, ρv, ρw)
    state.ρe = ρe_tot
    state.moisture.ρq_tot = FT(0)
    ρχ = zero(FT)
    if z <= 100
        ρχ += FT(0.1) * (cospi(z / 2 / 100))^2
    end
    state.tracers.ρχ = SVector{1, FT}(ρχ)
end

# Define problem configuration kernel
function config_problem(FT, N, resolution, xmax, ymax, zmax)

    ## Boundary conditions
    T_bot = FT(299)

    _cp_d::FT = cp_d(param_set)
    _grav::FT = grav(param_set)

    T_lapse = FT(_grav / _cp_d)
    T_top = T_bot - T_lapse * zmax

    ntracers = 1
    δ_χ = SVector{ntracers, FT}(1)

    ## Turbulence
    C_smag = FT(0.23)
    data_config = DryRayleighBenardConvectionDataConfig{FT}(
        0,
        0,
        0,
        xmax,
        ymax,
        zmax,
        T_bot,
        T_lapse,
        FT(T_bot - T_lapse * zmax),
    )

    ## Set up the model
    model = AtmosModel{FT}(
        AtmosLESConfigType,
        param_set;
        turbulence = Vreman(C_smag),
        source = (Gravity(),),
        boundarycondition = (
            AtmosBC(
                momentum = Impenetrable(NoSlip()),
                energy = PrescribedTemperature((state, aux, t) -> T_bot),
            ),
            AtmosBC(
                momentum = Impenetrable(NoSlip()),
                energy = PrescribedTemperature((state, aux, t) -> T_top),
            ),
        ),
        tracers = NTracers{ntracers, FT}(δ_χ),
        init_state_prognostic = init_problem!,
        data_config = data_config,
    )

    ## Set up the time-integrator, using a multirate infinitesimal step
    ## method. The option `splitting_type = ClimateMachine.SlowFastSplitting()`
    ## separates fast-slow modes by splitting away the acoustic waves and
    ## treating them via a sub-stepped explicit method.
    ode_solver = ClimateMachine.MISSolverType(;
        splitting_type = ClimateMachine.SlowFastSplitting(),
        mis_method = MIS2,
        fast_method = LSRK144NiegemannDiehlBusch,
        nsubsteps = 10,
    )

    config = ClimateMachine.AtmosLESConfiguration(
        "SO_cfSite_Test_Simulation",
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        param_set,
        init_problem!,
        solver_type = ode_solver,
        model = model,
    )
    return config
end

# Define diagnostics configuration kernel
function config_diagnostics(driver_config)
    interval = "10000steps"
    dgngrp = setup_atmos_default_diagnostics(
        AtmosLESConfigType(),
        interval,
        driver_config.name,
    )
    return ClimateMachine.DiagnosticsConfiguration([dgngrp])
end

# Define main entry point kernel
function main()
    FT = Float64
    ## DG polynomial order
    N = 4
    ## Domain resolution and size
    Δh = FT(10)                                                                  # this is resolution?? what units?
    ## Time integrator setup
    t0 = FT(0)
    ## Courant number
    CFLmax = FT(20)
    timeend = FT(1000)
    xmax, ymax, zmax = FT(250), FT(250), FT(500)

    @testset "SO_cfSite_Test_Simulation" begin
        for Δh in Δh
            Δv = Δh
            resolution = (Δh, Δh, Δv)
            driver_config = config_problem(FT, N, resolution, xmax, ymax, zmax)
            solver_config = ClimateMachine.SolverConfiguration(
                t0,
                timeend,
                driver_config,
                init_on_cpu = true,
                Courant_number = CFLmax,
            )
            dgn_config = config_diagnostics(driver_config)
            ## User defined callbacks (TMAR positivity preserving filter)
            cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
                Filters.apply!(
                    solver_config.Q,
                    ("moisture.ρq_tot",),
                    solver_config.dg.grid,
                    TMARFilter(),
                )
                nothing
            end
            result = ClimateMachine.invoke!(
                solver_config;
                diagnostics_config = dgn_config,
                user_callbacks = (cbtmarfilter,),
                check_euclidean_distance = true,
            )
            ## result == engf/eng0
            @test isapprox(result, FT(1); atol = 1.5e-2)
        end
    end
end

# Run
main()
