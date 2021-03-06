export Dirichlet, Neumann

boundary_conditions(::LandModel) = (1, 2)

function boundary_state!(
    nf,
    bctype,
    land::LandModel,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    t,
    _...,
)
    args = (state⁺, diff⁺, aux⁺, n̂, state⁻, diff⁻, aux⁻, t)
    soil_boundary_state!(nf, bctype, land, land.soil, land.soil.water, args...)
    soil_boundary_state!(nf, bctype, land, land.soil, land.soil.heat, args...)
end


function boundary_state!(
    nf,
    bctype,
    land::LandModel,
    state⁺::Vars,
    aux⁺::Vars,
    n̂,
    state⁻::Vars,
    aux⁻::Vars,
    t,
    _...,
)
    args = (state⁺, aux⁺, n̂, state⁻, aux⁻, t)
    soil_boundary_state!(nf, bctype, land, land.soil, land.soil.water, args...)
    soil_boundary_state!(nf, bctype, land, land.soil, land.soil.heat, args...)
end

abstract type AbstractBoundaryFunctions end

"""
    struct Dirichlet{Fs, Fb} <: AbstractBoundaryFunctions

A concrete type to hold the surface state and bottom state variable
values/functions, if Dirichlet boundary conditions are desired.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Dirichlet{Fs, Fb} <: AbstractBoundaryFunctions
    "Surface state boundary condition"
    surface_state::Fs = nothing
    "Bottom state boundary condition"
    bottom_state::Fb = nothing
end

"""
    struct Neumann{Fs, Fb} <: AbstractBoundaryFunctions

A concrete type to hold the surface and/or bottom diffusive flux
values/functions, if Neumann boundary conditions are desired.

Note that these are intended to be scalar values. In the boundary_state!
functions, they are multiplied by the `ẑ` vector (i.e. the normal vector `n̂`
to the domain at the upper boundary, and -`n̂` at the lower boundary. These
normal vectors point out of the domain.)

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Neumann{Fs, Fb} <: AbstractBoundaryFunctions
    "Surface flux boundary condition"
    surface_flux::Fs = nothing
    "Bottom flux boundary condition"
    bottom_flux::Fb = nothing
end
