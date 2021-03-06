using Test
using MPI
using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using LinearAlgebra
using Random

include("advection_diffusion_model.jl")

struct Box{dim} end
struct Sphere end

vertical_unit_vector(::Box{2}, ::SVector{3}) = SVector(0, 1, 0)
vertical_unit_vector(::Box{3}, ::SVector{3}) = SVector(0, 0, 1)
vertical_unit_vector(::Sphere, coord::SVector{3}) = coord / norm(coord)

projection(::EveryDirection, ::SVector{3}) = I
projection(::VerticalDirection, k::SVector{3}) = k * k'
projection(::HorizontalDirection, k::SVector{3}) = I - k * k'

struct TestProblem{adv, diff, dir, topo} <: AdvectionDiffusionProblem end

initial_ρ(::Box, x) = prod(sin.(π * x))
function initial_ρ(::Sphere, x)
    r = norm(x)
    φ = atan(x[2], x[1])
    θ = atan(sqrt(x[1]^2 + x[2]^2), x[3])
    return sin(π * (r - 1)) * sin(φ) * sin(θ)
end

velocity(::Box, x) = sin.(π * x)
function velocity(::Sphere, x)
    r = norm(x)
    φ = atan(x[2], x[1])
    θ = atan(sqrt(x[1]^2 + x[2]^2), x[3])
    return sin(π * (r - 1)) .*
           SVector(cos(φ) * cos(θ), sin(φ) * cos(θ), cos(φ) * sin(θ))
end

function init_velocity_diffusion!(
    ::TestProblem{adv, diff, dir, topo},
    aux::Vars,
    geom::LocalGeometry,
) where {adv, diff, dir, topo}
    k = vertical_unit_vector(topo, geom.coord)
    P = projection(dir, k)
    aux.advection.u =
        !isnothing(adv) ? P * velocity(topo, geom.coord) : zeros(SVector{3})
    aux.diffusion.D =
        !isnothing(diff) ? SMatrix{3, 3}(P) / 200 : zeros(SMatrix{3, 3})
end

function initial_condition!(
    ::TestProblem{adv, diff, dir, topo},
    state,
    aux,
    localgeo,
    t,
) where {adv, diff, dir, topo}
    state.ρ = initial_ρ(topo, localgeo.coord)
end

function create_topology(::Box{dim}, mpicomm, Ne, FT) where {dim}
    brickrange = ntuple(j -> range(FT(0); length = Ne + 1, stop = 1), dim)
    periodicity = ntuple(j -> false, dim)
    bc = ntuple(j -> (3, 3), dim)
    StackedBrickTopology(
        mpicomm,
        brickrange;
        periodicity = periodicity,
        boundary = bc,
    )
end

function create_topology(::Sphere, mpicomm, Ne, FT)
    vert_range = grid1d(FT(1), FT(2), nelem = Ne)
    StackedCubedSphereTopology(mpicomm, Ne, vert_range, boundary = (3, 3))
end

create_dg(model, grid, direction) = DGModel(
    model,
    grid,
    RusanovNumericalFlux(),
    CentralNumericalFluxSecondOrder(),
    CentralNumericalFluxGradient(),
    direction = direction,
)

function test_run(
    adv,
    diff,
    topo,
    mpicomm,
    ArrayType,
    FT,
    polynomialorder,
    Ne,
    level,
)
    topology = create_topology(topo(), mpicomm, Ne, FT)
    grid = DiscontinuousSpectralElementGrid(
        topology,
        FloatType = FT,
        DeviceArray = ArrayType,
        polynomialorder = polynomialorder,
        meshwarp = topo == Sphere ? cubedshellwarp : (x...) -> identity(x),
    )

    problems = (
        p = TestProblem{adv, diff, EveryDirection(), topo()}(),
        vp = TestProblem{adv, diff, VerticalDirection(), topo()}(),
        hp = TestProblem{adv, diff, HorizontalDirection(), topo()}(),
    )

    models = map(AdvectionDiffusion{3}, problems)

    dgmodels = map(models) do m
        (
            dg = create_dg(m, grid, EveryDirection()),
            vdg = create_dg(m, grid, VerticalDirection()),
            hdg = create_dg(m, grid, HorizontalDirection()),
        )
    end

    Q = init_ode_state(dgmodels.p.dg, FT(0), init_on_cpu = true)
    # do one Euler step to trigger numerical fluxes in subsequent evaluations
    let
        dt = 1e-3
        dQ = similar(Q)
        dgmodels.p.dg(dQ, Q, nothing, FT(0))
        Q .+= dt .* dQ
    end

    # evaluate all combinations
    dQ = map(
        x -> map(dg -> (dQ = similar(Q); dg(dQ, Q, nothing, FT(0)); dQ), x),
        dgmodels,
    )

    # set up tolerances
    atolm = 6e-13
    atolv = 3e-4 / 5^(level - 1)
    atolh = 0.0003 / 5^(level - 1)

    @testset "total" begin
        atol = topo <: Box || diff == nothing ? atolm : atolh
        @test isapprox(norm(dQ.p.dg .- dQ.p.vdg .- dQ.p.hdg), 0, atol = atol)
        @test isapprox(norm(dQ.vp.dg .- dQ.vp.vdg .- dQ.vp.hdg), 0, atol = atol)
        @test isapprox(norm(dQ.hp.dg .- dQ.hp.vdg .- dQ.hp.hdg), 0, atol = atol)
    end

    @testset "vertical" begin
        atol = topo <: Box ? atolm : atolv
        @test isapprox(norm(dQ.vp.dg .- dQ.vp.vdg), 0, atol = atol)
        @test isapprox(norm(dQ.vp.hdg), 0, atol = atol)
        @test isapprox(norm(dQ.vp.dg .- dQ.p.vdg), 0, atol = atol)
    end

    @testset "horizontal" begin
        atol = topo <: Box ? atolm : atolh
        @test isapprox(norm(dQ.hp.dg .- dQ.hp.hdg), 0, atol = atol)
        @test isapprox(norm(dQ.hp.vdg), 0, atol = atol)
        @test isapprox(norm(dQ.hp.dg - dQ.p.hdg), 0, atol = atol)
    end
end

let
    ClimateMachine.init()
    ArrayType = ClimateMachine.array_type()
    mpicomm = MPI.COMM_WORLD
    FT = Float64
    numlevels = 2
    polynomialorder = 4
    base_num_elem = 4

    @testset "$(@__FILE__)" begin
        @testset for topo in (Box{2}, Box{3}, Sphere)
            @testset for (adv, diff) in (
                (Advection, nothing),
                (nothing, Diffusion),
                (Advection, Diffusion),
            )
                @testset for level in 1:numlevels
                    Ne = 2^(level - 1) * base_num_elem
                    test_run(
                        adv,
                        diff,
                        topo,
                        mpicomm,
                        ArrayType,
                        FT,
                        polynomialorder,
                        Ne,
                        level,
                    )
                end
            end
        end
    end
end
nothing
