# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - STEP calculation

"""
    StepPoint

Result of equilibrium calculation at a single point in STEP calculation.

# Fields

  - `axis_value::Float64`: Value of swept variable (T or x)
  - `gibbs_energy::Float64`: Total Gibbs energy [J/mol]
  - `n_phases::Int`: Number of stable phases (1 or 2)
  - `phase_amounts::Vector{Float64}`: Amount of each phase [mol]
  - `phase_compositions::Vector{Vector{Float64}}`: Composition of each phase
  - `converged::Bool`: Whether calculation converged
"""
struct StepPoint
    axis_value::Float64
    gibbs_energy::Float64
    n_phases::Int
    phase_amounts::Vector{Float64}
    phase_compositions::Vector{Vector{Float64}}
    converged::Bool
end

"""
    StepResult

Complete result of a STEP calculation.

# Fields

  - `axis_variable::Symbol`: Variable being swept (:T or :X)
  - `axis_values::Vector{Float64}`: Values of swept variable
  - `points::Vector{StepPoint}`: Results at each point
  - `x_overall::Float64`: Fixed overall composition (first component mole fraction)
"""
struct StepResult
    axis_variable::Symbol
    axis_values::Vector{Float64}
    points::Vector{StepPoint}
    x_overall::Float64
end

"""
    step_temperature(phase::Phase, db::Database, x_overall, T_start, T_stop, T_step;
                     P=1e5, component_indices=(1,2))

Perform STEP calculation over temperature range.

# Arguments

  - `phase::Phase`: Phase to calculate
  - `db::Database`: Thermodynamic database
  - `x_overall::Real`: Overall composition (second component mole fraction)
  - `T_start::Real`: Start temperature [K]
  - `T_stop::Real`: Stop temperature [K]
  - `T_step::Real`: Temperature step [K]
  - `P::Real=1e5`: Pressure [Pa]
  - `component_indices::Tuple{Int,Int}=(1,2)`: Component indices in sublattice 1.
    For multi-component databases (e.g., steel1.TDB), specify the desired components:
    e.g., (1, 3) for Cr-Mo where Cr=1, Mo=3 in BCC_A2.

# Returns

  - `StepResult`: Results at each temperature point

# Example

```julia
# Standard binary (AgCu)
result = step_temperature(fcc, db, 0.7, 800.0, 1000.0, 100.0)

# Multi-component (Cr-Mo from steel1.TDB)
result = step_temperature(bcc, db, 0.5, 1800.0, 2000.0, 100.0; component_indices = (1, 3))
```
"""
function step_temperature(
    phase::Phase,
    db::Database,
    x_overall::Real,
    T_start::Real,
    T_stop::Real,
    T_step::Real;
    P::Real = 1e5,
    component_indices::Tuple{Int, Int} = (1, 2),
)
    # Generate temperature array
    T_values = collect(T_start:T_step:T_stop)

    # Calculate equilibrium at each temperature
    points = StepPoint[]
    for T in T_values
        point = calculate_step_point(phase, T, db, x_overall, P, component_indices)
        push!(points, point)
    end

    return StepResult(:T, T_values, points, Float64(x_overall))
end

"""
    calculate_step_point(phase, T, db, x_overall, P, component_indices)

Calculate equilibrium at a single temperature point.
"""
function calculate_step_point(
    phase::Phase,
    T::Real,
    db::Database,
    x_overall::Real,
    P::Real,
    component_indices::Tuple{Int, Int} = (1, 2),
)
    # Try two-phase equilibrium first
    # Note: calculate_miscibility_equilibrium will be updated to support component_indices
    # For now, it uses default (1,2), which may not work for multi-component databases
    result = calculate_miscibility_equilibrium(phase, T, db, x_overall; P = P)

    if !isnothing(result)
        # Two-phase: extract from TwoPhaseResult
        G_total =
            result.phase1_fraction * result.phase1.gibbs_energy +
            result.phase2_fraction * result.phase2.gibbs_energy

        # Extract compositions
        i1, _ = component_indices
        x1 = result.phase1.site_fractions[1, i1]  # First component in phase 1
        x2 = result.phase2.site_fractions[1, i1]  # First component in phase 2

        return StepPoint(
            Float64(T),
            G_total,
            2,
            [result.phase1_fraction, result.phase2_fraction],
            [[x1, 1 - x1], [x2, 1 - x2]],
            true,
        )
    else
        # Single-phase: calculate G directly at given composition using make_y_matrix
        y = make_y_matrix(phase, x_overall, component_indices)

        G = calculate_gibbs_energy(phase, T, y, db; P = P)

        return StepPoint(
            Float64(T),
            G,
            1,
            [1.0],
            [[1 - x_overall, x_overall]],  # [first_component, second_component]
            true,
        )
    end
end

# Helper functions for plotting

"""
    temperatures(result::StepResult) -> Vector{Float64}

Return array of temperature values from STEP result.
"""
temperatures(result::StepResult) = result.axis_values

"""
    gibbs_energies(result::StepResult) -> Vector{Float64}

Return array of Gibbs energy values from STEP result.
"""
gibbs_energies(result::StepResult) = [p.gibbs_energy for p in result.points]

"""
    phase_fractions(result::StepResult, phase_idx::Int) -> Vector{Float64}

Return array of phase fraction values for specified phase index.
"""
function phase_fractions(result::StepResult, phase_idx::Int)
    return [
        phase_idx <= p.n_phases ? p.phase_amounts[phase_idx] : 0.0
        for p in result.points
    ]
end

"""
    Base.show(io::IO, result::StepResult)

Pretty print STEP result.
"""
function Base.show(io::IO, result::StepResult)
    n = length(result.points)
    T_min = minimum(result.axis_values)
    T_max = maximum(result.axis_values)
    println(io, "StepResult: $(n) points")
    println(io, "  Variable: $(result.axis_variable)")
    println(io, "  Range: $(T_min) to $(T_max)")
    println(io, "  x_overall: $(result.x_overall)")
end

# Fallback methods for type errors

"""
    step_temperature(phase, db, x_overall, T_start, T_stop, T_step; P=1e5)

Fallback method for invalid phase type.
"""
function step_temperature(phase, db, x_overall, T_start, T_stop, T_step; kwargs...)
    error("Invalid phase type: $(typeof(phase)). Expected Phase.")
end

"""
    step_temperature(phase::Phase, db, x_overall, T_start, T_stop, T_step; P=1e5)

Fallback method for invalid database type.
"""
function step_temperature(phase::Phase, db, x_overall, T_start, T_stop, T_step; kwargs...)
    error("Invalid database type: $(typeof(db)). Expected Database.")
end
