# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - MAP calculation for phase diagrams

"""
    PhaseDiagramPoint

Single point on a phase diagram boundary.

# Fields

  - `temperature::Float64`: Temperature [K]
  - `compositions::NTuple{2, Float64}`: Boundary compositions (left, right) as mole fractions
  - `converged::Bool`: Whether the calculation converged

# Note

`compositions[1]` is the left (lower) boundary, `compositions[2]` is the right (upper) boundary.
"""
struct PhaseDiagramPoint
    temperature::Float64
    compositions::NTuple{2, Float64}
    converged::Bool
end

"""
    PhaseDiagramResult

Complete result of a MAP (phase diagram) calculation.

# Fields

  - `phase::Phase`: The phase being analyzed
  - `points::Vector{PhaseDiagramPoint}`: Boundary points at each temperature
"""
struct PhaseDiagramResult
    phase::Phase
    points::Vector{PhaseDiagramPoint}
end

"""
    map_phase_diagram(phase::Phase, db::Database, T_start, T_stop, T_step; P=1e5, component_indices=(1,2))

Calculate phase diagram (MAP) by finding phase boundaries over a temperature range.

For phases with miscibility gaps, this finds the two-phase region boundaries
at each temperature point.

# Arguments

  - `phase::Phase`: Phase to analyze
  - `db::Database`: Thermodynamic database
  - `T_start::Real`: Start temperature [K]
  - `T_stop::Real`: Stop temperature [K]
  - `T_step::Real`: Temperature step [K]
  - `P::Real=1e5`: Pressure [Pa]
  - `component_indices::Tuple{Int,Int}=(1,2)`: Indices of binary components in first sublattice

# Returns

  - `PhaseDiagramResult`: Phase boundaries at each temperature
"""
function map_phase_diagram(
    phase::Phase,
    db::Database,
    T_start::Real,
    T_stop::Real,
    T_step::Real;
    P::Real = 1e5,
    component_indices::Tuple{Int, Int} = (1, 2),
)
    # Generate temperature array
    T_values = collect(T_start:T_step:T_stop)

    # Calculate phase boundaries at each temperature
    points = PhaseDiagramPoint[]
    for T in T_values
        point = calculate_map_point(phase, T, db, P, component_indices)
        push!(points, point)
    end

    return PhaseDiagramResult(phase, points)
end

"""
    calculate_map_point(phase, T, db, P, component_indices)

Calculate phase boundaries at a single temperature.
"""
function calculate_map_point(
    phase::Phase,
    T::Real,
    db::Database,
    P::Real,
    component_indices::Tuple{Int, Int},
)
    # Find miscibility gap at this temperature
    gap = find_miscibility_gap(phase, T, db; P = P, component_indices = component_indices)

    if isnothing(gap)
        # No miscibility gap found - return unconverged point
        return PhaseDiagramPoint(Float64(T), (0.0, 0.0), false)
    end

    # Ensure compositions[1] < compositions[2]
    x1, x2 = minmax(gap.x1, gap.x2)

    return PhaseDiagramPoint(Float64(T), (x1, x2), true)
end

# Helper functions for plotting

"""
    temperatures(result::PhaseDiagramResult) -> Vector{Float64}

Return array of temperature values from MAP result.
"""
temperatures(result::PhaseDiagramResult) = [p.temperature for p in result.points]

"""
    boundaries(result::PhaseDiagramResult, idx::Int) -> Vector{Float64}

Return array of boundary compositions at index `idx`.

# Arguments

  - `result::PhaseDiagramResult`: MAP calculation result
  - `idx::Int`: Boundary index (1 = left/lower, 2 = right/upper)

# Examples

```julia
result = map_phase_diagram(fcc, db, 800.0, 1200.0, 50.0)
left_boundary = boundaries(result, 1)   # Left (lower) compositions
right_boundary = boundaries(result, 2)  # Right (upper) compositions
```
"""
function boundaries(result::PhaseDiagramResult, idx::Int)
    if idx < 1 || idx > 2
        throw(ArgumentError("idx must be 1 or 2, got $idx"))
    end
    return [p.compositions[idx] for p in result.points]
end

"""
    Base.show(io::IO, result::PhaseDiagramResult)

Pretty print MAP result.
"""
function Base.show(io::IO, result::PhaseDiagramResult)
    n = length(result.points)
    n_converged = count(p -> p.converged, result.points)
    T_min = minimum(p.temperature for p in result.points)
    T_max = maximum(p.temperature for p in result.points)
    println(io, "PhaseDiagramResult: $(n) points ($(n_converged) converged)")
    println(io, "  Phase: $(result.phase.name)")
    println(io, "  T range: $(T_min) to $(T_max) K")
end

# Fallback methods for type errors

"""
    map_phase_diagram(phase, db, T_start, T_stop, T_step; P=1e5)

Fallback method for invalid phase type.
"""
function map_phase_diagram(phase, db, T_start, T_stop, T_step; kwargs...)
    error("Invalid phase type: $(typeof(phase)). Expected Phase.")
end

"""
    map_phase_diagram(phase::Phase, db, T_start, T_stop, T_step; P=1e5)

Fallback method for invalid database type.
"""
function map_phase_diagram(phase::Phase, db, T_start, T_stop, T_step; kwargs...)
    error("Invalid database type: $(typeof(db)). Expected Database.")
end
