# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Two-phase equilibrium calculations

using ForwardDiff

"""
    TwoPhaseResult

Result of two-phase equilibrium calculation.

# Fields

  - `phase1::SinglePhaseResult`: First phase result
  - `phase2::SinglePhaseResult`: Second phase result
  - `phase1_fraction::Float64`: Molar fraction of phase 1
  - `phase2_fraction::Float64`: Molar fraction of phase 2
  - `chemical_potentials::Vector{Float64}`: Equilibrium chemical potentials
  - `converged::Bool`: Whether the calculation converged
  - `iterations::Int`: Number of iterations
"""
struct TwoPhaseResult
    phase1::SinglePhaseResult
    phase2::SinglePhaseResult
    phase1_fraction::Float64
    phase2_fraction::Float64
    chemical_potentials::Vector{Float64}
    converged::Bool
    iterations::Int
end

"""
    find_common_tangent(phase1::Phase, phase2::Phase, T::Real, db::Database;
                        P::Real=1e5, n_points::Int=101, tol::Float64=1e-6)
                        -> Union{TwoPhaseResult, Nothing}

Find two-phase equilibrium using common tangent construction.

For binary systems, finds compositions where chemical potentials are equal
in both phases.

# Arguments

  - `phase1`, `phase2`: The two phases
  - `T`: Temperature [K]
  - `db`: Database
  - `P`: Pressure [Pa]
  - `n_points`: Grid points for initial search
  - `tol`: Tolerance for chemical potential equality

# Returns

  - `TwoPhaseResult` if equilibrium found, `nothing` otherwise
"""
function find_common_tangent(
    phase1::Phase,
    phase2::Phase,
    T::Real,
    db::Database;
    P::Real = 1e5,
    n_points::Int = 101,
    tol::Float64 = 1e-6,
)
    # Grid scan both phases
    solver = GridSearchSolver(n_points = n_points)
    result1 = scan_composition(phase1, T, db, solver; P = P)
    result2 = scan_composition(phase2, T, db, solver; P = P)

    # Compute chemical potentials at each grid point
    diff_method = ForwardDiffMethod()
    mu1_all = Vector{Vector{Float64}}()
    mu2_all = Vector{Vector{Float64}}()

    for y in result1.y_values
        mu = compute_chemical_potential(phase1, T, y, db, diff_method; P = P)
        push!(mu1_all, mu)
    end

    for y in result2.y_values
        mu = compute_chemical_potential(phase2, T, y, db, diff_method; P = P)
        push!(mu2_all, mu)
    end

    # Find closest chemical potential match
    best_i, best_j = 1, 1
    best_diff = Inf

    for i in eachindex(mu1_all)
        for j in eachindex(mu2_all)
            # Check if both components have equal chemical potentials
            diff = sum(abs.(mu1_all[i] .- mu2_all[j]))
            if diff < best_diff
                best_diff = diff
                best_i, best_j = i, j
            end
        end
    end

    # Refine with Newton iteration if a good candidate was found
    if best_diff < 1e4  # Reasonable starting point
        y1_init = result1.y_values[best_i]
        y2_init = result2.y_values[best_j]

        result = refine_two_phase_equilibrium(
            phase1, phase2, T, db, y1_init, y2_init;
            P = P, tol = tol,
        )

        if !isnothing(result)
            return result
        end
    end

    return nothing
end

"""
    refine_two_phase_equilibrium(phase1::Phase, phase2::Phase, T::Real, db::Database,
                                  y1_init::Matrix{Float64}, y2_init::Matrix{Float64};
                                  P::Real=1e5, tol::Float64=1e-6, max_iter::Int=100)
                                  -> Union{TwoPhaseResult, Nothing}

Refine two-phase equilibrium using Newton iteration.

Solves the system:
μ₁^α - μ₁^β = 0
μ₂^α - μ₂^β = 0 (automatically satisfied if first is satisfied for binary)

For binary with site fraction y, the independent variable is y₁ for each phase.
"""
function refine_two_phase_equilibrium(
    phase1::Phase,
    phase2::Phase,
    T::Real,
    db::Database,
    y1_init::Matrix{Float64},
    y2_init::Matrix{Float64};
    P::Real = 1e5,
    tol::Float64 = 1e-6,
    max_iter::Int = 100,
)
    diff_method = ForwardDiffMethod()

    # For binary single-sublattice: y1[1,1] and y2[1,1] are the independent variables
    y1 = copy(y1_init)
    y2 = copy(y2_init)

    min_frac = 1e-8

    converged = false
    iterations = 0

    for iter in 1:max_iter
        iterations = iter

        # Compute chemical potentials
        mu1 = compute_chemical_potential(phase1, T, y1, db, diff_method; P = P)
        mu2 = compute_chemical_potential(phase2, T, y2, db, diff_method; P = P)

        # Residual: difference in chemical potentials
        # For binary: μ₁^α - μ₁^β = 0 and μ₂^α - μ₂^β = 0
        # But μ₁ - μ₂ is related to dG/dx, so we use Δμ = μ₁^α - μ₂^α - (μ₁^β - μ₂^β)
        delta_mu1 = mu1[1] - mu1[2]  # Chemical potential difference in phase 1
        delta_mu2 = mu2[1] - mu2[2]  # Chemical potential difference in phase 2

        res = delta_mu1 - delta_mu2

        if abs(res) < tol
            converged = true
            break
        end

        # Compute Jacobian numerically
        h = 1e-6

        # d(res)/d(y1[1,1])
        y1_p = copy(y1)
        y1_p[1, 1] = min(1 - min_frac, y1[1, 1] + h)
        y1_p[1, 2] = 1 - y1_p[1, 1]
        mu1_p = compute_chemical_potential(phase1, T, y1_p, db, diff_method; P = P)
        dres_dy1 = ((mu1_p[1] - mu1_p[2]) - delta_mu1) / h

        # d(res)/d(y2[1,1])
        y2_p = copy(y2)
        y2_p[1, 1] = min(1 - min_frac, y2[1, 1] + h)
        y2_p[1, 2] = 1 - y2_p[1, 1]
        mu2_p = compute_chemical_potential(phase2, T, y2_p, db, diff_method; P = P)
        dres_dy2 = -((mu2_p[1] - mu2_p[2]) - delta_mu2) / h

        # We have 1 equation, 2 unknowns
        # Add constraint: adjust both equally to find tangent
        # Newton step: minimize |res|^2 subject to moving along tangent line
        # Simple approach: move both y1 and y2 to reduce residual

        J_norm = sqrt(dres_dy1^2 + dres_dy2^2)
        if J_norm < 1e-15
            break
        end

        # Gradient descent direction
        alpha = 0.5
        step = -res / J_norm

        dy1 = alpha * step * dres_dy1 / J_norm
        dy2 = alpha * step * dres_dy2 / J_norm

        # Update with bounds
        y1[1, 1] = clamp(y1[1, 1] + dy1, min_frac, 1 - min_frac)
        y1[1, 2] = 1 - y1[1, 1]
        y2[1, 1] = clamp(y2[1, 1] + dy2, min_frac, 1 - min_frac)
        y2[1, 2] = 1 - y2[1, 1]
    end

    if !converged
        return nothing
    end

    # Compute final results
    G1 = calculate_gibbs_energy(phase1, T, y1, db; P = P)
    G2 = calculate_gibbs_energy(phase2, T, y2, db; P = P)
    mu_final = compute_chemical_potential(phase1, T, y1, db, diff_method; P = P)

    result1 = SinglePhaseResult(phase1, y1, G1, true, iterations)
    result2 = SinglePhaseResult(phase2, y2, G2, true, iterations)

    # Phase fractions are undetermined without overall composition
    # Return 0.5/0.5 as default
    return TwoPhaseResult(result1, result2, 0.5, 0.5, mu_final, true, iterations)
end

"""
    calculate_two_phase_equilibrium(phase1::Phase, phase2::Phase, T::Real, db::Database,
                                     x_overall::Float64;
                                     P::Real=1e5, solver::EquilibriumSolver=HybridSolver())
                                     -> Union{TwoPhaseResult, Nothing}

Calculate two-phase equilibrium at given overall composition.

# Arguments

  - `phase1`, `phase2`: The two phases
  - `T`: Temperature [K]
  - `db`: Database
  - `x_overall`: Overall mole fraction of first component
  - `P`: Pressure [Pa]
  - `solver`: Solver for single-phase optimization

# Returns

  - `TwoPhaseResult` if two phases are stable, `nothing` if single phase
"""
function calculate_two_phase_equilibrium(
    phase1::Phase,
    phase2::Phase,
    T::Real,
    db::Database,
    x_overall::Float64;
    P::Real = 1e5,
    solver::EquilibriumSolver = HybridSolver(),
)
    # First find the common tangent
    ct_result = find_common_tangent(phase1, phase2, T, db; P = P)

    if isnothing(ct_result)
        return nothing
    end

    # Get equilibrium compositions
    x1 = ct_result.phase1.site_fractions[1, 1]  # Composition in phase 1
    x2 = ct_result.phase2.site_fractions[1, 1]  # Composition in phase 2

    # Check if x_overall is within the two-phase region
    x_min, x_max = minmax(x1, x2)

    if x_overall < x_min || x_overall > x_max
        # Single phase region
        return nothing
    end

    # Lever rule to determine phase fractions
    if abs(x1 - x2) < 1e-10
        # Same composition - no two-phase region
        return nothing
    end

    f1 = (x_overall - x2) / (x1 - x2)
    f2 = 1 - f1

    # Clamp to valid range
    f1 = clamp(f1, 0.0, 1.0)
    f2 = 1 - f1

    return TwoPhaseResult(
        ct_result.phase1,
        ct_result.phase2,
        f1,
        f2,
        ct_result.chemical_potentials,
        true,
        ct_result.iterations,
    )
end

"""
    find_phase_boundaries(phase1::Phase, phase2::Phase, T_range::AbstractRange,
                          db::Database; P::Real=1e5)
                          -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}

Find phase boundaries over a temperature range.

# Returns

  - `temperatures`: Temperature points where equilibrium was found
  - `x1_boundary`: Composition boundary of phase 1
  - `x2_boundary`: Composition boundary of phase 2
"""
function find_phase_boundaries(
    phase1::Phase,
    phase2::Phase,
    T_range::AbstractRange,
    db::Database;
    P::Real = 1e5,
)
    temperatures = Float64[]
    x1_boundary = Float64[]
    x2_boundary = Float64[]

    for T in T_range
        result = find_common_tangent(phase1, phase2, T, db; P = P, n_points = 51)

        if !isnothing(result) && result.converged
            push!(temperatures, T)
            push!(x1_boundary, result.phase1.site_fractions[1, 1])
            push!(x2_boundary, result.phase2.site_fractions[1, 1])
        end
    end

    return temperatures, x1_boundary, x2_boundary
end

"""
    MiscibilityGapResult

Result of miscibility gap calculation within a single phase.

# Fields

  - `phase::Phase`: The phase with miscibility gap
  - `x1::Float64`: Composition of first region (mole fraction of first component)
  - `x2::Float64`: Composition of second region
  - `G1::Float64`: Gibbs energy at x1
  - `G2::Float64`: Gibbs energy at x2
  - `converged::Bool`: Whether the calculation converged
  - `iterations::Int`: Number of iterations
"""
struct MiscibilityGapResult
    phase::Phase
    x1::Float64
    x2::Float64
    G1::Float64
    G2::Float64
    converged::Bool
    iterations::Int
end

"""
    find_miscibility_gap(phase::Phase, T::Real, db::Database,
                         x1_init::Float64=0.05, x2_init::Float64=0.95;
                         P::Real=1e5, max_iter::Int=100, tol::Float64=1e-8,
                         component_indices::Tuple{Int,Int}=(1,2))
                         -> Union{MiscibilityGapResult, Nothing}

Find miscibility gap (spinodal decomposition) within a single phase using
common tangent construction with Newton-Raphson refinement.

For binary systems, finds two compositions x1 and x2 where the common tangent
condition is satisfied: dG/dx|_{x1} = dG/dx|_{x2} = (G2-G1)/(x2-x1)

# Arguments

  - `phase`: Phase to analyze
  - `T`: Temperature [K]
  - `db`: Database
  - `x1_init`: Initial guess for left composition (default: 0.05)
  - `x2_init`: Initial guess for right composition (default: 0.95)
  - `P`: Pressure [Pa]
  - `max_iter`: Maximum Newton iterations
  - `tol`: Convergence tolerance
  - `component_indices`: Tuple of (first, second) component indices in sublattice 1.
    Default is (1, 2). For multi-component databases, specify the desired components.

# Returns

  - `MiscibilityGapResult` if miscibility gap found, `nothing` otherwise

# Example

```julia
# Standard binary (AgCu)
gap = find_miscibility_gap(fcc, 1000.0, db)

# Multi-component (Cr-Mo from steel1.TDB, CR=2, MO=4 in LIQUID)
gap = find_miscibility_gap(liquid, 2000.0, db; component_indices = (2, 4))
```
"""
function find_miscibility_gap(
    phase::Phase,
    T::Real,
    db::Database,
    x1_init::Float64 = 0.05,
    x2_init::Float64 = 0.95;
    P::Real = 1e5,
    max_iter::Int = 100,
    tol::Float64 = 1e-8,
    component_indices::Tuple{Int, Int} = (1, 2),
)
    i1, i2 = component_indices

    # Helper function: G as function of x (first component mole fraction)
    # Note: x represents the mole fraction of component at i1 (first component)
    function G_at_x(x)
        n_const = maximum(length.(phase.constituents))
        n_sub = length(phase.sites)
        # Use promoted type to support ForwardDiff Dual numbers
        y = fill(zero(x), n_sub, n_const)
        y[1, i1] = x            # First component
        y[1, i2] = one(x) - x   # Second component
        # Fill other sublattices with VA
        for s in 2:n_sub
            va_idx = find_va_index(phase.constituents[s])
            y[s, va_idx] = one(x)
        end
        return calculate_gibbs_energy(phase, T, y, db; P = P)
    end

    function dG_at_x(x)
        DI.derivative(G_at_x, DI.AutoForwardDiff(), x)
    end

    function d2G_at_x(x)
        DI.second_derivative(G_at_x, DI.AutoForwardDiff(), x)
    end

    x1, x2 = x1_init, x2_init
    converged = false
    iterations = 0

    for iter in 1:max_iter
        iterations = iter

        G1, G2 = G_at_x(x1), G_at_x(x2)
        dG1, dG2 = dG_at_x(x1), dG_at_x(x2)
        d2G1, d2G2 = d2G_at_x(x1), d2G_at_x(x2)

        # Common tangent slope
        Δx = x2 - x1
        if abs(Δx) < 1e-10
            break
        end
        slope = (G2 - G1) / Δx

        # Residuals: dG/dx should equal slope at both points
        r1 = dG1 - slope
        r2 = dG2 - slope

        err = sqrt(r1^2 + r2^2)
        if err < tol
            converged = true
            break
        end

        # Jacobian matrix for Newton-Raphson
        J11 = d2G1 + slope / Δx - dG1 / Δx
        J12 = -dG2 / Δx + slope / Δx
        J21 = dG1 / Δx - slope / Δx
        J22 = d2G2 - slope / Δx + dG2 / Δx

        # Solve J * Δ = -r
        det = J11 * J22 - J12 * J21
        if abs(det) < 1e-15
            break
        end

        Δx1 = (-J22 * r1 + J12 * r2) / det
        Δx2 = (J21 * r1 - J11 * r2) / det

        # Line search with backtracking
        alpha = 1.0
        for _ in 1:10
            x1_new = clamp(x1 + alpha * Δx1, 0.001, 0.499)
            x2_new = clamp(x2 + alpha * Δx2, 0.501, 0.999)

            G1_new, G2_new = G_at_x(x1_new), G_at_x(x2_new)
            dG1_new, dG2_new = dG_at_x(x1_new), dG_at_x(x2_new)
            slope_new = (G2_new - G1_new) / (x2_new - x1_new)

            r1_new = dG1_new - slope_new
            r2_new = dG2_new - slope_new
            err_new = sqrt(r1_new^2 + r2_new^2)

            if err_new < err
                x1, x2 = x1_new, x2_new
                break
            end
            alpha *= 0.5
        end
    end

    if !converged
        return nothing
    end

    G1, G2 = G_at_x(x1), G_at_x(x2)
    return MiscibilityGapResult(phase, x1, x2, G1, G2, true, iterations)
end

"""
    find_miscibility_gap(phase::Phase, T::Real, x1_init, x2_init; kwargs...) -> Union{MiscibilityGapResult, Nothing}

Find miscibility gap using Julia DSL parameters (no Database required).

This method is used when the phase has `julia_parameters` defined.
It enables parameter optimization with automatic differentiation.

# Arguments

  - `phase`: Phase with julia_parameters defined
  - `T`: Temperature [K]
  - `x1_init`: Initial guess for composition 1 (default: 0.05)
  - `x2_init`: Initial guess for composition 2 (default: 0.95)
  - `P`: Pressure [Pa] (keyword)
  - `max_iter`: Maximum Newton iterations (keyword)
  - `tol`: Convergence tolerance (keyword)

# Returns

  - `MiscibilityGapResult` if miscibility gap found, `nothing` otherwise
"""
function find_miscibility_gap(
    phase::Phase,
    T::Real,
    x1_init::Float64 = 0.05,
    x2_init::Float64 = 0.95;
    P::Real = 1e5,
    max_iter::Int = 100,
    tol::Float64 = 1e-8,
)
    if !has_julia_parameters(phase)
        error(
            "Phase $(phase.name) has no Julia parameters. Use find_miscibility_gap(phase, T, db) with a Database.",
        )
    end

    # Helper function: G as function of x (first component mole fraction)
    function G_at_x(x)
        n_const = length(phase.constituents[1])
        n_sub = length(phase.sites)
        # Use promoted type to support ForwardDiff Dual numbers
        y = fill(zero(x), n_sub, n_const)
        y[1, 1] = x
        y[1, 2] = one(x) - x
        # Fill other sublattices (typically vacancy)
        for s in 2:n_sub
            y[s, 1] = one(x)
        end
        return calculate_gibbs_energy(phase, T, y; P = P)
    end

    function dG_at_x(x)
        DI.derivative(G_at_x, DI.AutoForwardDiff(), x)
    end

    function d2G_at_x(x)
        DI.second_derivative(G_at_x, DI.AutoForwardDiff(), x)
    end

    x1, x2 = x1_init, x2_init
    converged = false
    iterations = 0

    for iter in 1:max_iter
        iterations = iter

        G1, G2 = G_at_x(x1), G_at_x(x2)
        dG1, dG2 = dG_at_x(x1), dG_at_x(x2)
        d2G1, d2G2 = d2G_at_x(x1), d2G_at_x(x2)

        # Common tangent slope
        Δx = x2 - x1
        if abs(Δx) < 1e-10
            break
        end
        slope = (G2 - G1) / Δx

        # Residuals: dG/dx should equal slope at both points
        r1 = dG1 - slope
        r2 = dG2 - slope

        err = sqrt(r1^2 + r2^2)
        if err < tol
            converged = true
            break
        end

        # Jacobian matrix for Newton-Raphson
        J11 = d2G1 + slope / Δx - dG1 / Δx
        J12 = -dG2 / Δx + slope / Δx
        J21 = dG1 / Δx - slope / Δx
        J22 = d2G2 - slope / Δx + dG2 / Δx

        # Solve J * Δ = -r
        det = J11 * J22 - J12 * J21
        if abs(det) < 1e-15
            break
        end

        Δx1 = (-J22 * r1 + J12 * r2) / det
        Δx2 = (J21 * r1 - J11 * r2) / det

        # Line search with backtracking
        alpha = 1.0
        for _ in 1:10
            x1_new = clamp(x1 + alpha * Δx1, 0.001, 0.499)
            x2_new = clamp(x2 + alpha * Δx2, 0.501, 0.999)

            G1_new, G2_new = G_at_x(x1_new), G_at_x(x2_new)
            dG1_new, dG2_new = dG_at_x(x1_new), dG_at_x(x2_new)
            slope_new = (G2_new - G1_new) / (x2_new - x1_new)

            r1_new = dG1_new - slope_new
            r2_new = dG2_new - slope_new
            err_new = sqrt(r1_new^2 + r2_new^2)

            if err_new < err
                x1, x2 = x1_new, x2_new
                break
            end
            alpha *= 0.5
        end
    end

    if !converged
        return nothing
    end

    G1, G2 = G_at_x(x1), G_at_x(x2)
    return MiscibilityGapResult(phase, x1, x2, G1, G2, true, iterations)
end

"""
    calculate_miscibility_equilibrium(phase::Phase, T::Real, db::Database,
                                       x_overall::Float64;
                                       P::Real=1e5, component_indices::Tuple{Int,Int}=(1,2),
                                       kwargs...)
                                       -> Union{TwoPhaseResult, Nothing}

Calculate equilibrium for a phase with miscibility gap at given overall composition.

# Arguments

  - `phase`: Phase with potential miscibility gap
  - `T`: Temperature [K]
  - `db`: Database
  - `x_overall`: Overall mole fraction of second component
  - `P`: Pressure [Pa]
  - `component_indices`: Tuple of (first, second) component indices in sublattice 1.

# Returns

  - `TwoPhaseResult` if in two-phase region, `nothing` if single phase
"""
function calculate_miscibility_equilibrium(
    phase::Phase,
    T::Real,
    db::Database,
    x_overall::Float64;
    P::Real = 1e5,
    component_indices::Tuple{Int, Int} = (1, 2),
    kwargs...,
)
    gap = find_miscibility_gap(
        phase,
        T,
        db;
        P = P,
        component_indices = component_indices,
        kwargs...,
    )

    if isnothing(gap)
        return nothing
    end

    # Check if x_overall is within the miscibility gap
    x_min, x_max = minmax(gap.x1, gap.x2)

    if x_overall < x_min || x_overall > x_max
        return nothing
    end

    # Lever rule
    f1 = (gap.x2 - x_overall) / (gap.x2 - gap.x1)
    f2 = 1 - f1

    # Create site fraction matrices
    # Note: gap.x1 and gap.x2 are first component mole fractions
    i1, i2 = component_indices
    n_const = maximum(length.(phase.constituents))
    n_sub = length(phase.sites)

    y1 = zeros(n_sub, n_const)
    y1[1, i1] = gap.x1
    y1[1, i2] = 1 - gap.x1
    for s in 2:n_sub
        va_idx = find_va_index(phase.constituents[s])
        y1[s, va_idx] = 1.0
    end

    y2 = zeros(n_sub, n_const)
    y2[1, i1] = gap.x2
    y2[1, i2] = 1 - gap.x2
    for s in 2:n_sub
        va_idx = find_va_index(phase.constituents[s])
        y2[s, va_idx] = 1.0
    end

    result1 = SinglePhaseResult(phase, y1, gap.G1, true, gap.iterations)
    result2 = SinglePhaseResult(phase, y2, gap.G2, true, gap.iterations)

    # Chemical potential at equilibrium
    diff_method = ForwardDiffMethod()
    mu = compute_chemical_potential(phase, T, y1, db, diff_method; P = P)

    return TwoPhaseResult(result1, result2, f1, f2, mu, true, gap.iterations)
end
