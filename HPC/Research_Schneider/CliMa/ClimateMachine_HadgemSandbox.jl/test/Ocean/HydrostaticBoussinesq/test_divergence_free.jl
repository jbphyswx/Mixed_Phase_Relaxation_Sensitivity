using Test
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Grids: polynomialorder
using ClimateMachine.BalanceLaws: vars_state, Prognostic
using ClimateMachine.Ocean.HydrostaticBoussinesq

using CLIMAParameters
using CLIMAParameters.Planet: grav
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

function config_simple_box(FT, N, resolution, dimensions; BC = nothing)
    if BC == nothing
        problem = HomogeneousBox{FT}(dimensions...)
    else
        problem = HomogeneousBox{FT}(dimensions...; BC = BC)
    end

    _grav::FT = grav(param_set)
    cʰ = sqrt(_grav * problem.H) # m/s
    model = HydrostaticBoussinesqModel{FT}(param_set, problem, cʰ = cʰ)

    config = ClimateMachine.OceanBoxGCMConfiguration(
        "homogeneous_box",
        N,
        resolution,
        param_set,
        model,
    )

    return config
end

function test_divergence_free(; imex::Bool = false, BC = nothing)
    FT = Float64

    # DG polynomial order
    N = Int(4)

    # Domain resolution and size
    Nˣ = Int(4)
    Nʸ = Int(4)
    Nᶻ = Int(4)
    resolution = (Nˣ, Nʸ, Nᶻ)

    Lˣ = 4e6   # m
    Lʸ = 4e6   # m
    H = 400   # m
    dimensions = (Lˣ, Lʸ, H)

    timestart = FT(0)    # s
    timeend = FT(36000) # s

    if imex
        solver_type = ClimateMachine.IMEXSolverType(
            implicit_model = LinearHBModel,
            implicit_solver = ClimateMachine.SystemSolvers.SingleColumnLU,
        )
        Courant_number = 0.1
    else
        solver_type = ClimateMachine.ExplicitSolverType(
            solver_method = LSRK144NiegemannDiehlBusch,
        )
        Courant_number = 0.4
    end

    driver_config = config_simple_box(FT, N, resolution, dimensions; BC = BC)

    grid = driver_config.grid
    vert_filter = CutoffFilter(grid, polynomialorder(grid) - 1)
    exp_filter = ExponentialFilter(grid, 1, 8)
    modeldata = (vert_filter = vert_filter, exp_filter = exp_filter)

    solver_config = ClimateMachine.SolverConfiguration(
        timestart,
        timeend,
        driver_config,
        init_on_cpu = true,
        ode_solver_type = solver_type,
        ode_dt = 600,
        modeldata = modeldata,
        Courant_number = Courant_number,
    )

    result = ClimateMachine.invoke!(solver_config)

    maxQ = Vars{vars_state(driver_config.bl, Prognostic(), FT)}(maximum(
        solver_config.Q,
        dims = (1, 3),
    ))
    minQ = Vars{vars_state(driver_config.bl, Prognostic(), FT)}(minimum(
        solver_config.Q,
        dims = (1, 3),
    ))

    @test maxQ.θ ≈ minQ.θ
end

@testset "$(@__FILE__)" begin
    boundary_conditions = [(
        OceanBC(Impenetrable(NoSlip()), Insulating()),
        OceanBC(Impenetrable(NoSlip()), Insulating()),
        OceanBC(Penetrable(KinematicStress()), Insulating()),
    )]

    for BC in boundary_conditions
        test_divergence_free(imex = false, BC = BC)
        test_divergence_free(imex = true, BC = BC)
    end
end
