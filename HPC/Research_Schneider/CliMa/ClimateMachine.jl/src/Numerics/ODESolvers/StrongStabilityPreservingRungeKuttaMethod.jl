export StrongStabilityPreservingRungeKutta
export SSPRK22Heuns, SSPRK22Ralstons, SSPRK33ShuOsher, SSPRK34SpiteriRuuth

"""
    StrongStabilityPreservingRungeKutta(f, RKA, RKB, RKC, Q; dt, t0 = 0)

This is a time stepping object for explicitly time stepping the differential
equation given by the right-hand-side function `f` with the state `Q`, i.e.,

```math
  \\dot{Q} = f(Q, t)
```

with the required time step size `dt` and optional initial time `t0`.  This
time stepping object is intended to be passed to the `solve!` command.

The constructor builds a strong-stability-preserving Runge--Kutta scheme
based on the provided `RKA`, `RKB` and `RKC` coefficient arrays.

The available concrete implementations are:

  - [`SSPRK33ShuOsher`](@ref)
  - [`SSPRK34SpiteriRuuth`](@ref)
"""
mutable struct StrongStabilityPreservingRungeKutta{T, RT, AT, Nstages} <:
               AbstractODESolver
    "time step"
    dt::RT
    "time"
    t::RT
    "elapsed time steps"
    steps::Int
    "rhs function"
    rhs!::Any
    "Storage for RHS during the `StrongStabilityPreservingRungeKutta` update"
    Rstage::AT
    "Storage for the stage state during the `StrongStabilityPreservingRungeKutta` update"
    Qstage::AT
    "RK coefficient vector A (rhs scaling)"
    RKA::Array{RT, 2}
    "RK coefficient vector B (rhs add in scaling)"
    RKB::Array{RT, 1}
    "RK coefficient vector C (time scaling)"
    RKC::Array{RT, 1}

    function StrongStabilityPreservingRungeKutta(
        rhs!,
        RKA,
        RKB,
        RKC,
        Q::AT;
        dt = 0,
        t0 = 0,
    ) where {AT <: AbstractArray}
        T = eltype(Q)
        RT = real(T)
        new{T, RT, AT, length(RKB)}(
            RT(dt),
            RT(t0),
            0,
            rhs!,
            similar(Q),
            similar(Q),
            RKA,
            RKB,
            RKC,
        )
    end
end

"""
    dostep!(Q, ssp::StrongStabilityPreservingRungeKutta, p, time::Real,
            nsteps::Int, iStage::Int, [slow_δ, slow_rv_dQ, slow_scaling])

Wrapper function to use the strong stability preserving Runge--Kutta method `ssp`
as the fast solver for a Multirate Infinitesimal Step method by calling dostep!(Q,
ssp::StrongStabilityPreservingRungeKutta, p, time::Real, [slow_δ, slow_rv_dQ,
slow_scaling]) nsubsteps times.
"""
function dostep!(
    Q,
    ssp::StrongStabilityPreservingRungeKutta,
    p,
    time::Real,
    nsteps::Int,
    iStage::Int,
    slow_δ = nothing,
    slow_rv_dQ = nothing,
    slow_scaling = nothing,
)
    for i in 1:nsteps
        dostep!(Q, ssp, p, time, slow_δ, slow_rv_dQ, slow_scaling)
        time += ssp.dt
    end
end

"""
    ODESolvers.dostep!(Q, ssp::StrongStabilityPreservingRungeKutta, p,
                       time::Real, [slow_δ, slow_rv_dQ, slow_scaling])

Use the strong stability preserving Runge--Kutta method `ssp` to step `Q`
forward in time from the current time `time` to final time `time + getdt(ssp)`.

If the optional parameter `slow_δ !== nothing` then `slow_rv_dQ * slow_δ` is
added as an additional ODE right-hand side source. If the optional parameter
`slow_scaling !== nothing` then after the final stage update the scaling
`slow_rv_dQ *= slow_scaling` is performed.
"""
function dostep!(
    Q,
    ssp::StrongStabilityPreservingRungeKutta,
    p,
    time,
    slow_δ = nothing,
    slow_rv_dQ = nothing,
    in_slow_scaling = nothing,
)
    dt = ssp.dt

    RKA, RKB, RKC = ssp.RKA, ssp.RKB, ssp.RKC
    rhs! = ssp.rhs!
    Rstage, Qstage = ssp.Rstage, ssp.Qstage

    rv_Q = realview(Q)
    rv_Rstage = realview(Rstage)
    rv_Qstage = realview(Qstage)
    groupsize = 256

    rv_Qstage .= rv_Q
    for s in 1:length(RKB)
        rhs!(Rstage, Qstage, p, time + RKC[s] * dt, increment = false)

        slow_scaling = nothing
        if s == length(RKB)
            slow_scaling = in_slow_scaling
        end
        event = Event(array_device(Q))
        event = update!(array_device(Q), groupsize)(
            rv_Rstage,
            rv_Q,
            rv_Qstage,
            RKA[s, 1],
            RKA[s, 2],
            RKB[s],
            dt,
            slow_δ,
            slow_rv_dQ,
            slow_scaling;
            ndrange = length(rv_Q),
            dependencies = (event,),
        )
        wait(array_device(Q), event)
    end
    rv_Q .= rv_Qstage
end

@kernel function update!(
    dQ,
    Q,
    Qstage,
    rka1,
    rka2,
    rkb,
    dt,
    slow_δ,
    slow_dQ,
    slow_scaling,
)
    i = @index(Global, Linear)
    @inbounds begin
        if slow_δ !== nothing
            dQ[i] += slow_δ * slow_dQ[i]
        end
        Qstage[i] = rka1 * Q[i] + rka2 * Qstage[i] + dt * rkb * dQ[i]
        if slow_scaling !== nothing
            slow_dQ[i] *= slow_scaling
        end
    end
end

"""
    SSPRK22Heuns(f, Q; dt, t0 = 0)

This function returns a [`StrongStabilityPreservingRungeKutta`](@ref) time stepping object
for explicitly time stepping the differential
equation given by the right-hand-side function `f` with the state `Q`, i.e.,

```math
  \\dot{Q} = f(Q, t)
```

with the required time step size `dt` and optional initial time `t0`.  This
time stepping object is intended to be passed to the `solve!` command.

This uses the second-order, 2-stage, strong-stability-preserving, Runge--Kutta scheme
of Shu and Osher (1988) (also known as Heun's method.)
Exact choice of coefficients from wikipedia page for Heun's method :)

### References
 - [Shu1988](@cite)
 - [Heun1900](@cite)
"""
function SSPRK22Heuns(F, Q::AT; dt = 0, t0 = 0) where {AT <: AbstractArray}
    T = eltype(Q)
    RT = real(T)
    RKA = [RT(1) RT(0); RT(1 // 2) RT(1 // 2)]
    RKB = [RT(1), RT(1 // 2)]
    RKC = [RT(0), RT(1)]
    StrongStabilityPreservingRungeKutta(F, RKA, RKB, RKC, Q; dt = dt, t0 = t0)
end

"""
    SSPRK22Ralstons(f, Q; dt, t0 = 0)

This function returns a [`StrongStabilityPreservingRungeKutta`](@ref) time stepping object
for explicitly time stepping the differential
equation given by the right-hand-side function `f` with the state `Q`, i.e.,

```math
  \\dot{Q} = f(Q, t)
```

with the required time step size `dt` and optional initial time `t0`.  This
time stepping object is intended to be passed to the `solve!` command.

This uses the second-order, 2-stage, strong-stability-preserving, Runge--Kutta scheme
of Shu and Osher (1988) (also known as Ralstons's method.)
Exact choice of coefficients from wikipedia page for Heun's method :)

### References
 - [Shu1988](@cite)
 - [Ralston1962](@cite)
"""
function SSPRK22Ralstons(F, Q::AT; dt = 0, t0 = 0) where {AT <: AbstractArray}
    T = eltype(Q)
    RT = real(T)
    RKA = [RT(1) RT(0); RT(5 // 8) RT(3 // 8)]
    RKB = [RT(2 // 3), RT(3 // 4)]
    RKC = [RT(0), RT(2 // 3)]
    StrongStabilityPreservingRungeKutta(F, RKA, RKB, RKC, Q; dt = dt, t0 = t0)
end

"""
    SSPRK33ShuOsher(f, Q; dt, t0 = 0)

This function returns a [`StrongStabilityPreservingRungeKutta`](@ref) time stepping object
for explicitly time stepping the differential
equation given by the right-hand-side function `f` with the state `Q`, i.e.,

```math
  \\dot{Q} = f(Q, t)
```

with the required time step size `dt` and optional initial time `t0`.  This
time stepping object is intended to be passed to the `solve!` command.

This uses the third-order, 3-stage, strong-stability-preserving, Runge--Kutta scheme
of Shu and Osher (1988)

### References
 - [Shu1988](@cite)
"""
function SSPRK33ShuOsher(F, Q::AT; dt = 0, t0 = 0) where {AT <: AbstractArray}
    T = eltype(Q)
    RT = real(T)
    RKA = [RT(1) RT(0); RT(3 // 4) RT(1 // 4); RT(1 // 3) RT(2 // 3)]
    RKB = [RT(1), RT(1 // 4), RT(2 // 3)]
    RKC = [RT(0), RT(1), RT(1 // 2)]
    StrongStabilityPreservingRungeKutta(F, RKA, RKB, RKC, Q; dt = dt, t0 = t0)
end

"""
    SSPRK34SpiteriRuuth(f, Q; dt, t0 = 0)

This function returns a [`StrongStabilityPreservingRungeKutta`](@ref) time stepping object
for explicitly time stepping the differential
equation given by the right-hand-side function `f` with the state `Q`, i.e.,

```math
  \\dot{Q} = f(Q, t)
```

with the required time step size `dt` and optional initial time `t0`.  This
time stepping object is intended to be passed to the `solve!` command.

This uses the third-order, 4-stage, strong-stability-preserving, Runge--Kutta scheme
of Spiteri and Ruuth (1988)

### References
 - [Spiteri2002](@cite)
"""
function SSPRK34SpiteriRuuth(
    F,
    Q::AT;
    dt = 0,
    t0 = 0,
) where {AT <: AbstractArray}
    T = eltype(Q)
    RT = real(T)
    RKA = [RT(1) RT(0); RT(0) RT(1); RT(2 // 3) RT(1 // 3); RT(0) RT(1)]
    RKB = [RT(1 // 2); RT(1 // 2); RT(1 // 6); RT(1 // 2)]
    RKC = [RT(0); RT(1 // 2); RT(1); RT(1 // 2)]
    StrongStabilityPreservingRungeKutta(F, RKA, RKB, RKC, Q; dt = dt, t0 = t0)
end
