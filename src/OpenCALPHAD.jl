# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl

module OpenCALPHAD

using LinearAlgebra
using Printf
using ForwardDiff
using StaticArrays
using DifferentiationInterface
import DifferentiationInterface as DI

# Constants
const R = 8.314462618  # Gas constant [J/(mol·K)]
const P_REF = 1.0e5    # Reference pressure [Pa]

# Upstream version tracking (OpenCALPHAD Fortran)
# Based on local clone at ../opencalphad
const UPSTREAM_REPO = "https://github.com/sundmanbo/opencalphad"
const UPSTREAM_COMMIT = "e857325f5d38e32e41513186ddf1d8ae4b2e0349"
const UPSTREAM_DATE = "2026-01-01"

"""
    upstream_version()

Return information about the OpenCALPHAD (Fortran) version this package is based on.

# Example

```julia
julia> upstream_version()
(repo = "https://github.com/sundmanbo/opencalphad", commit = "0e7c26e6...", date = "2026-01-23")
```
"""
function upstream_version()
    return (
        repo = UPSTREAM_REPO,
        commit = UPSTREAM_COMMIT,
        date = UPSTREAM_DATE,
    )
end

# Type definitions
include("Types/elements.jl")
include("Types/functions.jl")
include("Types/phases.jl")
include("Types/conditions.jl")
include("Types/equilibrium.jl")

# Database I/O
include("Database/database.jl")
include("Database/tdb_parser.jl")
include("Database/expression.jl")

# Thermodynamic models
include("Models/gibbs.jl")
include("Models/thermodynamics.jl")

# Minimizer - interface and implementations
include("Minimizer/interface.jl")
include("Minimizer/differentiation.jl")
include("Minimizer/grid_minimizer.jl")
include("Minimizer/newton_solver.jl")

# Calculations
include("Calculations/equilibrium.jl")
include("Calculations/two_phase.jl")
include("Calculations/step.jl")
include("Calculations/map.jl")
include("Calculations/phase_field.jl")

# Plotting
include("Plotting/recipes.jl")

# Exports - Types
export Element, Phase, Parameter, GFunction, Database
export JuliaParameter  # DSL support
export EquilCondition, ConditionType, ConditionSet, Equilibrium
export PhaseResult, ComponentResult

# Exports - Condition types
export COND_T, COND_P, COND_N, COND_X, COND_W, COND_MU, COND_AC, COND_NP, COND_BP, COND_VP

# Exports - Solver types
export DiffMethod, ForwardDiffMethod, FiniteDiffMethod, DIMethod
export EquilibriumSolver, GridSearchSolver, NewtonSolver, HybridSolver
export SinglePhaseResult, GridScanResult

# Exports - Functions
export read_tdb
export calculate_equilibrium, step_calculation
export calculate_gibbs_energy
export calculate_entropy, calculate_enthalpy
export calculate_mixing_gibbs, calculate_mixing_enthalpy, calculate_mixing_entropy
export thermodynamic_properties
export num_sublattices, total_sites, stable_phases, n_stable_phases
export set_condition!, degrees_of_freedom
export get_element, get_phase, get_function
export evaluate_gfunction, compile_gfunction!
export evaluate_parameter

# Exports - Minimizer functions
export minimize_gibbs
export scan_composition, scan_composition_multi
export find_minimum_envelope, find_two_phase_region
export compute_gradient, compute_hessian, compute_chemical_potential
export clamp_site_fractions!

# Exports - Two-phase equilibrium
export TwoPhaseResult, MiscibilityGapResult
export find_common_tangent, calculate_two_phase_equilibrium, find_phase_boundaries
export find_miscibility_gap, calculate_miscibility_equilibrium

# Exports - STEP calculation
export StepPoint, StepResult
export step_temperature
export temperatures, gibbs_energies, phase_fractions

# Exports - MAP calculation (phase diagram)
export PhaseDiagramPoint, PhaseDiagramResult
export map_phase_diagram
export boundaries

# Exports - Plotting
export plot_phase_diagram, plot_gibbs_curve

# Exports - DSL (Julia function parameters)
export set_G!, set_L!
export add_element!, add_elements!, add_phase!
export has_julia_parameters, is_endmember, is_interaction
export get_julia_endmembers, get_julia_interactions
export validate_ad_compatibility

# Exports - Phase Field integration
export driving_force, chemical_potential, diffusion_potential
export phase_field_params

# Exports - Version info
export upstream_version

end # module
