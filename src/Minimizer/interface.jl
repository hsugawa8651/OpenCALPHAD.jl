# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Abstract interfaces for equilibrium solvers

# Differentiation Methods

"""
    DiffMethod

Abstract type for differentiation methods.
"""
abstract type DiffMethod end

"""
    ForwardDiffMethod <: DiffMethod

Use ForwardDiff.jl for automatic differentiation.
Fast and accurate, but may fail with some expression types.
"""
struct ForwardDiffMethod <: DiffMethod end

"""
    FiniteDiffMethod <: DiffMethod

Use finite differences for numerical differentiation.
Always works, but slower and less accurate.

# Fields

  - `epsilon::Float64`: Step size for finite differences (default: 1e-8)
"""
struct FiniteDiffMethod <: DiffMethod
    epsilon::Float64
    FiniteDiffMethod(eps::Float64 = 1e-8) = new(eps)
end

"""
    DIMethod{B} <: DiffMethod

Use DifferentiationInterface.jl for automatic differentiation.
Supports multiple backends (ForwardDiff, Enzyme, etc.).

# Type Parameters

  - `B`: Backend type (e.g., `AutoForwardDiff`, `AutoEnzyme`)

# Fields

  - `backend::B`: The differentiation backend

# Example

```julia
method = DIMethod()                      # Default: ForwardDiff backend
method = DIMethod(DI.AutoForwardDiff())  # Explicit ForwardDiff
method = DIMethod(DI.AutoEnzyme())       # Enzyme backend (future)
```
"""
struct DIMethod{B} <: DiffMethod
    backend::B
end

# Default constructor uses ForwardDiff backend
DIMethod() = DIMethod(DI.AutoForwardDiff())#==============================================================================##==============================================================================#

# Equilibrium Solvers

"""
    EquilibriumSolver

Abstract type for equilibrium calculation solvers.
"""
abstract type EquilibriumSolver end

"""
    GridSearchSolver <: EquilibriumSolver

Grid search over composition space for initial guess.

# Fields

  - `n_points::Int`: Number of grid points per dimension (default: 21)
  - `min_fraction::Float64`: Minimum site fraction to avoid log(0) (default: 1e-8)
"""
struct GridSearchSolver <: EquilibriumSolver
    n_points::Int
    min_fraction::Float64

    function GridSearchSolver(; n_points::Int = 21, min_fraction::Float64 = 1e-8)
        new(n_points, min_fraction)
    end
end

"""
    NewtonSolver <: EquilibriumSolver

Newton-Raphson solver using NLsolve.jl.

# Fields

  - `diff_method::DiffMethod`: Differentiation method
  - `max_iter::Int`: Maximum iterations (default: 100)
  - `ftol::Float64`: Function tolerance (default: 1e-8)
  - `xtol::Float64`: Variable tolerance (default: 1e-8)
"""
struct NewtonSolver <: EquilibriumSolver
    diff_method::DiffMethod
    max_iter::Int
    ftol::Float64
    xtol::Float64

    function NewtonSolver(;
        diff_method::DiffMethod = ForwardDiffMethod(),
        max_iter::Int = 100,
        ftol::Float64 = 1e-8,
        xtol::Float64 = 1e-8,
    )
        new(diff_method, max_iter, ftol, xtol)
    end
end

"""
    HybridSolver <: EquilibriumSolver

Hybrid solver: Grid search for initial guess, then Newton refinement.

# Fields

  - `grid::GridSearchSolver`: Grid search solver for initial guess
  - `newton::NewtonSolver`: Newton solver for refinement
"""
struct HybridSolver <: EquilibriumSolver
    grid::GridSearchSolver
    newton::NewtonSolver

    function HybridSolver(;
        n_points::Int = 11,
        diff_method::DiffMethod = ForwardDiffMethod(),
        max_iter::Int = 100,
        ftol::Float64 = 1e-8,
    )
        grid = GridSearchSolver(n_points = n_points)
        newton = NewtonSolver(diff_method = diff_method, max_iter = max_iter, ftol = ftol)
        new(grid, newton)
    end
end#==============================================================================##==============================================================================#

# Result Types

"""
    SinglePhaseResult

Result of single-phase Gibbs energy minimization.

# Fields

  - `phase::Phase`: The phase
  - `site_fractions::Matrix{Float64}`: Optimized site fractions y[sublattice, constituent]
  - `gibbs_energy::Float64`: Gibbs energy at minimum [J/mol]
  - `converged::Bool`: Whether the solver converged
  - `iterations::Int`: Number of iterations
"""
struct SinglePhaseResult
    phase::Phase
    site_fractions::Matrix{Float64}
    gibbs_energy::Float64
    converged::Bool
    iterations::Int
end

"""
    GridScanResult

Result of grid scan over composition space.

# Fields

  - `phase::Phase`: The phase
  - `x_grid::Vector{Float64}`: Composition grid points
  - `G_values::Vector{Float64}`: Gibbs energy at each grid point
  - `y_values::Vector{Matrix{Float64}}`: Site fractions at each grid point
  - `min_idx::Int`: Index of minimum G point
"""
struct GridScanResult
    phase::Phase
    x_grid::Vector{Float64}
    G_values::Vector{Float64}
    y_values::Vector{Matrix{Float64}}
    min_idx::Int
end#==============================================================================##==============================================================================#

# Generic Interface Functions (to be implemented by specific solvers)

"""
    minimize_gibbs(phase::Phase, T::Real, db::Database, solver::EquilibriumSolver;
                   x_target::Union{Nothing, Dict}=nothing, P::Real=1e5)

Minimize Gibbs energy for a single phase.

# Arguments

  - `phase`: Phase to minimize
  - `T`: Temperature [K]
  - `db`: Database with parameters
  - `solver`: Solver to use
  - `x_target`: Target composition (optional, for constrained minimization)
  - `P`: Pressure [Pa]

# Returns

  - `SinglePhaseResult`: Optimization result
"""
function minimize_gibbs end

"""
    scan_composition(phase::Phase, T::Real, db::Database, solver::GridSearchSolver;
                     element1::String="", element2::String="", P::Real=1e5)

Scan Gibbs energy over composition space for a binary system.

# Arguments

  - `phase`: Phase to scan
  - `T`: Temperature [K]
  - `db`: Database
  - `solver`: Grid search solver
  - `element1`, `element2`: Elements defining the composition axis
  - `P`: Pressure [Pa]

# Returns

  - `GridScanResult`: Scan result
"""
function scan_composition end

"""
    compute_gradient(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                     method::DiffMethod; P::Real=1e5)

Compute gradient of Gibbs energy with respect to site fractions.

# Returns

  - `Vector{Float64}`: ∂G/∂y (flattened)
"""
function compute_gradient end

"""
    compute_hessian(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                    method::DiffMethod; P::Real=1e5)

Compute Hessian of Gibbs energy with respect to site fractions.

# Returns

  - `Matrix{Float64}`: ∂²G/∂y²
"""
function compute_hessian end#==============================================================================##==============================================================================#

# Fallback methods for unsupported types

"""
Fallback for unsupported solver types.
"""
function minimize_gibbs(phase::Phase, T::Real, db::Database, solver; kwargs...)
    error("Unsupported solver type: $(typeof(solver))")
end

"""
Fallback for unsupported differentiation methods.
"""
function compute_gradient(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
    method; P::Real = 1e5)
    error("Unsupported differentiation method: $(typeof(method))")
end

function compute_hessian(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
    method; P::Real = 1e5)
    error("Unsupported differentiation method: $(typeof(method))")
end
