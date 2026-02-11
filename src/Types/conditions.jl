# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Equilibrium condition definitions

"""
    ConditionType

Types of equilibrium conditions.
"""
@enum ConditionType begin
    COND_T          # Temperature [K]
    COND_P          # Pressure [Pa]
    COND_N          # Total moles
    COND_X          # Mole fraction
    COND_W          # Mass fraction
    COND_MU         # Chemical potential
    COND_AC         # Activity
    COND_NP         # Moles of phase
    COND_BP         # Mass of phase
    COND_VP         # Volume fraction of phase
end

"""
    EquilCondition

Represents a single equilibrium condition.
(Named EquilCondition to avoid conflict with Base.Condition)

# Fields

  - `type::ConditionType`: Type of condition
  - `component::String`: Component name (for X, W, MU, AC) or phase name (for NP, BP, VP)
  - `value::Float64`: Condition value
"""
struct EquilCondition
    type::ConditionType
    component::String
    value::Float64
end

# Convenience constructors
EquilCondition(type::ConditionType, value::Real) = EquilCondition(type, "", Float64(value))

function Base.show(io::IO, c::EquilCondition)
    if c.component == ""
        print(io, "$(c.type)=$(c.value)")
    else
        print(io, "$(c.type)($(c.component))=$(c.value)")
    end
end

"""
    ConditionSet

Collection of equilibrium conditions.
"""
mutable struct ConditionSet
    temperature::Float64      # T [K]
    pressure::Float64         # P [Pa]
    total_moles::Float64      # N [mol]
    compositions::Dict{String, Float64}  # X(component) or W(component)
    fixed_potentials::Dict{String, Float64}  # MU(component)
    fixed_phases::Dict{String, Float64}      # NP(phase)
end

function ConditionSet()
    return ConditionSet(
        298.15,     # Default T
        1.0e5,      # Default P
        1.0,        # Default N
        Dict{String, Float64}(),
        Dict{String, Float64}(),
        Dict{String, Float64}(),
    )
end

function Base.show(io::IO, cs::ConditionSet)
    print(io, "Conditions: T=$(cs.temperature)K, P=$(cs.pressure)Pa, N=$(cs.total_moles)")
    for (comp, x) in cs.compositions
        print(io, ", X($comp)=$x")
    end
end

"""
    set_condition!(cs::ConditionSet, cond::EquilCondition)

Add or update a condition in the condition set.
"""
function set_condition!(cs::ConditionSet, cond::EquilCondition)
    if cond.type == COND_T
        cs.temperature = cond.value
    elseif cond.type == COND_P
        cs.pressure = cond.value
    elseif cond.type == COND_N
        cs.total_moles = cond.value
    elseif cond.type == COND_X
        cs.compositions[cond.component] = cond.value
    elseif cond.type == COND_MU
        cs.fixed_potentials[cond.component] = cond.value
    elseif cond.type == COND_NP
        cs.fixed_phases[cond.component] = cond.value
    end
    return cs
end

"""
    degrees_of_freedom(cs::ConditionSet, n_components::Int) -> Int

Calculate the degrees of freedom for the system.
DOF = n_components + 2 - n_conditions
"""
function degrees_of_freedom(cs::ConditionSet, n_components::Int)
    n_conds = 2  # T and P are always set
    n_conds += length(cs.compositions)
    n_conds += length(cs.fixed_potentials)
    n_conds += length(cs.fixed_phases)
    return n_components + 2 - n_conds
end
