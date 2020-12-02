using Test, MPI
include("../testhelpers.jl")

@testset "ODE Solvers" begin
    tests = [(1, "ode_tests_basic.jl")]
    runmpi(tests, @__FILE__)
end
