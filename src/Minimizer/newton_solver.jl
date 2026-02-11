# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Newton-Raphson solver for Gibbs energy minimization

"""
    minimize_gibbs(phase::Phase, T::Real, db::Database, solver::NewtonSolver;
                   x_target::Union{Nothing, Dict}=nothing, P::Real=1e5) -> SinglePhaseResult

Minimize Gibbs energy using Newton-Raphson method with projected gradient.

For sublattice model with constraint Σy_s = 1 for each sublattice, we project
the gradient and Hessian onto the constraint manifold.
"""
function minimize_gibbs(
    phase::Phase,
    T::Real,
    db::Database,
    solver::NewtonSolver;
    x_target::Union{Nothing, Dict} = nothing,
    P::Real = 1e5,
)
    n_sub = length(phase.sites)
    n_const = length(phase.constituents[1])

    # Initial guess: equal fractions
    y = ones(n_sub, n_const) ./ n_const
    clamp_site_fractions!(y; min_val = solver.ftol)

    converged = false
    iterations = 0
    G = calculate_gibbs_energy(phase, T, y, db; P = P)

    for iter in 1:solver.max_iter
        iterations = iter

        # Compute gradient and Hessian
        grad = compute_gradient(phase, T, y, db, solver.diff_method; P = P)
        hess = compute_hessian(phase, T, y, db, solver.diff_method; P = P)

        # For single sublattice binary: reduce to 1D problem
        # y1 is the only free variable (y2 = 1 - y1)
        if n_sub == 1 && n_const == 2
            # Effective gradient: dG/dy1 considering y2 = 1 - y1
            # dG_eff = ∂G/∂y1 - ∂G/∂y2
            g_eff = grad[1] - grad[2]

            # Effective Hessian: d²G/dy1² considering constraint
            # H_eff = ∂²G/∂y1² - 2*∂²G/∂y1∂y2 + ∂²G/∂y2²
            H_eff = hess[1, 1] - 2 * hess[1, 2] + hess[2, 2]

            # Newton step
            if abs(H_eff) > 1e-15
                dy1 = -g_eff / H_eff
            else
                dy1 = -sign(g_eff) * 0.1  # Gradient descent fallback
            end

            # Line search with backtracking
            alpha = 1.0
            y_new = copy(y)
            G_new = G

            for _ in 1:20
                y_new[1, 1] = y[1, 1] + alpha * dy1
                y_new[1, 2] = 1.0 - y_new[1, 1]

                # Check bounds
                if y_new[1, 1] < solver.ftol
                    y_new[1, 1] = solver.ftol
                    y_new[1, 2] = 1.0 - y_new[1, 1]
                elseif y_new[1, 1] > 1.0 - solver.ftol
                    y_new[1, 1] = 1.0 - solver.ftol
                    y_new[1, 2] = 1.0 - y_new[1, 1]
                end

                G_new = calculate_gibbs_energy(phase, T, y_new, db; P = P)

                # Armijo condition
                if G_new <= G - 0.0001 * alpha * abs(g_eff * dy1)
                    break
                end
                alpha *= 0.5
            end

            # Update
            y .= y_new
            G_old = G
            G = G_new

            # Check convergence
            if abs(g_eff) < solver.ftol || abs(G - G_old) < solver.ftol
                converged = true
                break
            end
        else
            # General case: projected Newton for multiple sublattices
            # Use simple gradient descent with projection for now
            step_size = 0.01

            for s in 1:n_sub
                n_c = length(phase.constituents[s])
                idx_start = (s - 1) * n_const + 1
                idx_end = idx_start + n_c - 1

                grad_s = grad[idx_start:idx_end]

                # Project gradient onto constraint plane (sum = 0)
                grad_s_proj = grad_s .- _mean(grad_s)

                # Update
                y[s, 1:n_c] .-= step_size * grad_s_proj
            end

            clamp_site_fractions!(y; min_val = solver.ftol)
            G_old = G
            G = calculate_gibbs_energy(phase, T, y, db; P = P)

            # Check convergence
            if norm(grad) < solver.ftol * length(grad) || abs(G - G_old) < solver.ftol
                converged = true
                break
            end
        end
    end

    return SinglePhaseResult(phase, y, G, converged, iterations)
end

"""
    minimize_gibbs(phase::Phase, T::Real, db::Database, solver::HybridSolver;
                   x_target::Union{Nothing, Dict}=nothing, P::Real=1e5) -> SinglePhaseResult

Hybrid solver: Grid search for initial guess, then Newton refinement.
"""
function minimize_gibbs(
    phase::Phase,
    T::Real,
    db::Database,
    solver::HybridSolver;
    x_target::Union{Nothing, Dict} = nothing,
    P::Real = 1e5,
)
    # Step 1: Grid search for initial guess
    grid_result = scan_composition(phase, T, db, solver.grid; P = P)
    y_init = grid_result.y_values[grid_result.min_idx]

    n_sub = length(phase.sites)
    n_const = length(phase.constituents[1])

    # Initialize Newton solver with grid result
    y = copy(y_init)
    clamp_site_fractions!(y; min_val = solver.newton.ftol)

    converged = false
    iterations = 0
    G = calculate_gibbs_energy(phase, T, y, db; P = P)

    for iter in 1:solver.newton.max_iter
        iterations = iter

        grad = compute_gradient(phase, T, y, db, solver.newton.diff_method; P = P)
        hess = compute_hessian(phase, T, y, db, solver.newton.diff_method; P = P)

        if n_sub == 1 && n_const == 2
            g_eff = grad[1] - grad[2]
            H_eff = hess[1, 1] - 2 * hess[1, 2] + hess[2, 2]

            if abs(H_eff) > 1e-15
                dy1 = -g_eff / H_eff
            else
                dy1 = -sign(g_eff) * 0.1
            end

            alpha = 1.0
            y_new = copy(y)
            G_new = G

            for _ in 1:20
                y_new[1, 1] = y[1, 1] + alpha * dy1
                y_new[1, 2] = 1.0 - y_new[1, 1]

                if y_new[1, 1] < solver.newton.ftol
                    y_new[1, 1] = solver.newton.ftol
                    y_new[1, 2] = 1.0 - y_new[1, 1]
                elseif y_new[1, 1] > 1.0 - solver.newton.ftol
                    y_new[1, 1] = 1.0 - solver.newton.ftol
                    y_new[1, 2] = 1.0 - y_new[1, 1]
                end

                G_new = calculate_gibbs_energy(phase, T, y_new, db; P = P)

                if G_new <= G - 0.0001 * alpha * abs(g_eff * dy1)
                    break
                end
                alpha *= 0.5
            end

            y .= y_new
            G_old = G
            G = G_new

            if abs(g_eff) < solver.newton.ftol || abs(G - G_old) < solver.newton.ftol
                converged = true
                break
            end
        else
            step_size = 0.01
            for s in 1:n_sub
                n_c = length(phase.constituents[s])
                idx_start = (s - 1) * n_const + 1
                idx_end = idx_start + n_c - 1
                grad_s = grad[idx_start:idx_end]
                grad_s_proj = grad_s .- _mean(grad_s)
                y[s, 1:n_c] .-= step_size * grad_s_proj
            end

            clamp_site_fractions!(y; min_val = solver.newton.ftol)
            G_old = G
            G = calculate_gibbs_energy(phase, T, y, db; P = P)

            if norm(grad) < solver.newton.ftol * length(grad) ||
               abs(G - G_old) < solver.newton.ftol
                converged = true
                break
            end
        end
    end

    return SinglePhaseResult(phase, y, G, converged, iterations)
end

# Use LinearAlgebra for norm
using LinearAlgebra: norm

# Simple mean function to avoid Statistics dependency
_mean(x) = sum(x) / length(x)
