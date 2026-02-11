# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Equilibrium calculation module

"""
    calculate_equilibrium(db::Database, conditions::ConditionSet) -> Equilibrium

Calculate thermodynamic equilibrium for the given database and conditions.

# Algorithm

 1. Extract conditions (T, P, composition)
 2. Global grid minimization for initial phase estimate
 3. Newton-Raphson iteration for precise solution
 4. Post-processing and result formatting

# Arguments

  - `db::Database`: Thermodynamic database
  - `conditions::ConditionSet`: Equilibrium conditions

# Returns

  - `Equilibrium`: Calculated equilibrium state
"""
function calculate_equilibrium(db::Database, conditions::ConditionSet)
    eq = Equilibrium(conditions)

    T = conditions.temperature
    P = conditions.pressure

    # Validate conditions
    if degrees_of_freedom(conditions, length(db.elements) - 1) != 0
        @warn "Degrees of freedom ≠ 0, calculation may not converge"
    end

    # TODO: Implement full equilibrium calculation
    # 1. Global grid minimization
    # 2. Newton-Raphson refinement
    # 3. Phase stability tests

    eq.converged = false
    eq.iterations = 0

    return eq
end

"""
    calculate_equilibrium(db::Database;
                          T::Real=298.15,
                          P::Real=1e5,
                          N::Real=1.0,
                          x::Dict{String,Float64}=Dict()) -> Equilibrium

Convenience function for equilibrium calculation with keyword arguments.
"""
function calculate_equilibrium(db::Database;
    T::Real = 298.15,
    P::Real = 1e5,
    N::Real = 1.0,
    x::Dict{String, Float64} = Dict{String, Float64}())
    conds = ConditionSet()
    conds.temperature = Float64(T)
    conds.pressure = Float64(P)
    conds.total_moles = Float64(N)

    for (comp, val) in x
        conds.compositions[uppercase(comp)] = Float64(val)
    end

    return calculate_equilibrium(db, conds)
end

"""
    step_calculation(db::Database, conditions::ConditionSet,
                     variable::Symbol, range::AbstractRange) -> Vector{Equilibrium}

Perform stepping calculation over a range of values for a variable.

# Arguments

  - `db::Database`: Thermodynamic database
  - `conditions::ConditionSet`: Base conditions
  - `variable::Symbol`: Variable to step (:T, :P, :X, etc.)
  - `range`: Range of values to step through

# Returns

  - `Vector{Equilibrium}`: Equilibrium results for each step
"""
function step_calculation(db::Database, conditions::ConditionSet,
    variable::Symbol, range::AbstractRange)
    results = Equilibrium[]

    for val in range
        conds = deepcopy(conditions)

        if variable == :T
            conds.temperature = Float64(val)
        elseif variable == :P
            conds.pressure = Float64(val)
        end

        eq = calculate_equilibrium(db, conds)
        push!(results, eq)
    end

    return results
end
