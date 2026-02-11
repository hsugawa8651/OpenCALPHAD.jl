# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Differentiation methods for Gibbs energy calculations

"""
    compute_gradient(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                     method::ForwardDiffMethod; P::Real=1e5) -> Vector{Float64}

Compute gradient using ForwardDiff automatic differentiation.
"""
function compute_gradient(
    phase::Phase,
    T::Real,
    y::AbstractMatrix,
    db::Database,
    method::ForwardDiffMethod;
    P::Real = 1e5,
)
    y_flat = vec(y)
    y_size = size(y)

    grad = ForwardDiff.gradient(y_flat) do yv
        y_mat = reshape(yv, y_size)
        calculate_gibbs_energy(phase, T, y_mat, db; P = P)
    end

    return grad
end

"""
    compute_gradient(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                     method::FiniteDiffMethod; P::Real=1e5) -> Vector{Float64}

Compute gradient using finite differences.
"""
function compute_gradient(
    phase::Phase,
    T::Real,
    y::AbstractMatrix,
    db::Database,
    method::FiniteDiffMethod;
    P::Real = 1e5,
)
    y_flat = vec(y)
    y_size = size(y)
    n = length(y_flat)
    h = method.epsilon

    G0 = calculate_gibbs_energy(phase, T, y, db; P = P)
    grad = zeros(n)

    for i in 1:n
        y_p = copy(y_flat)
        y_p[i] += h
        G_p = calculate_gibbs_energy(phase, T, reshape(y_p, y_size), db; P = P)
        grad[i] = (G_p - G0) / h
    end

    return grad
end

"""
    compute_gradient(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                     method::DIMethod; P::Real=1e5) -> Vector{Float64}

Compute gradient using DifferentiationInterface.jl.
"""
function compute_gradient(
    phase::Phase,
    T::Real,
    y::AbstractMatrix,
    db::Database,
    method::DIMethod;
    P::Real = 1e5,
)
    y_flat = vec(y)
    y_size = size(y)

    G_func(yv) = calculate_gibbs_energy(phase, T, reshape(yv, y_size), db; P = P)
    grad = DI.gradient(G_func, method.backend, y_flat)

    return grad
end

"""
    compute_hessian(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                    method::ForwardDiffMethod; P::Real=1e5) -> Matrix{Float64}

Compute Hessian using ForwardDiff automatic differentiation.
"""
function compute_hessian(
    phase::Phase,
    T::Real,
    y::AbstractMatrix,
    db::Database,
    method::ForwardDiffMethod;
    P::Real = 1e5,
)
    y_flat = vec(y)
    y_size = size(y)

    hess = ForwardDiff.hessian(y_flat) do yv
        y_mat = reshape(yv, y_size)
        calculate_gibbs_energy(phase, T, y_mat, db; P = P)
    end

    return hess
end

"""
    compute_hessian(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                    method::FiniteDiffMethod; P::Real=1e5) -> Matrix{Float64}

Compute Hessian using finite differences.
"""
function compute_hessian(
    phase::Phase,
    T::Real,
    y::AbstractMatrix,
    db::Database,
    method::FiniteDiffMethod;
    P::Real = 1e5,
)
    y_flat = vec(y)
    y_size = size(y)
    n = length(y_flat)
    h = method.epsilon

    hess = zeros(n, n)

    for i in 1:n
        for j in i:n
            y_pp = copy(y_flat)
            y_pm = copy(y_flat)
            y_mp = copy(y_flat)
            y_mm = copy(y_flat)

            y_pp[i] += h
            y_pp[j] += h
            y_pm[i] += h
            y_pm[j] -= h
            y_mp[i] -= h
            y_mp[j] += h
            y_mm[i] -= h
            y_mm[j] -= h

            G_pp = calculate_gibbs_energy(phase, T, reshape(y_pp, y_size), db; P = P)
            G_pm = calculate_gibbs_energy(phase, T, reshape(y_pm, y_size), db; P = P)
            G_mp = calculate_gibbs_energy(phase, T, reshape(y_mp, y_size), db; P = P)
            G_mm = calculate_gibbs_energy(phase, T, reshape(y_mm, y_size), db; P = P)

            hess[i, j] = (G_pp - G_pm - G_mp + G_mm) / (4 * h^2)
            hess[j, i] = hess[i, j]
        end
    end

    return hess
end

"""
    compute_hessian(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                    method::DIMethod; P::Real=1e5) -> Matrix{Float64}

Compute Hessian using DifferentiationInterface.jl.
"""
function compute_hessian(
    phase::Phase,
    T::Real,
    y::AbstractMatrix,
    db::Database,
    method::DIMethod;
    P::Real = 1e5,
)
    y_flat = vec(y)
    y_size = size(y)

    G_func(yv) = calculate_gibbs_energy(phase, T, reshape(yv, y_size), db; P = P)
    hess = DI.hessian(G_func, method.backend, y_flat)

    return hess
end

"""
    compute_chemical_potential(phase::Phase, T::Real, y::AbstractMatrix, db::Database,
                               method::DiffMethod; P::Real=1e5) -> Vector{Float64}

Compute chemical potentials from Gibbs energy derivatives.

For sublattice model:
μᵢ = G + Σⱼ (δᵢⱼ - xⱼ) · ∂G/∂xⱼ

Returns chemical potential for each component [J/mol].
"""
function compute_chemical_potential(
    phase::Phase,
    T::Real,
    y::AbstractMatrix,
    db::Database,
    method::DiffMethod;
    P::Real = 1e5,
)
    G = calculate_gibbs_energy(phase, T, y, db; P = P)
    grad = compute_gradient(phase, T, y, db, method; P = P)

    # For simple substitutional solution (single sublattice):
    # μᵢ = G + ∂G/∂yᵢ - Σⱼ yⱼ · ∂G/∂yⱼ
    #    = G + ∂G/∂yᵢ - y · grad

    n_sub = size(y, 1)
    n_const = size(y, 2)
    y_flat = vec(y)

    # Simple case: assume first sublattice contains the components
    mu = zeros(n_const)
    grad_sum = dot(y_flat, grad)

    for i in 1:n_const
        idx = i  # Index in flattened array for sublattice 1
        mu[i] = G + grad[idx] - grad_sum
    end

    return mu
end

"""
    clamp_site_fractions!(y::AbstractMatrix; min_val::Float64=1e-10, normalize::Bool=true)

Clamp site fractions to avoid numerical issues with log(0).
Optionally renormalize each sublattice to sum to 1.
Iterates to ensure min_val constraint holds after normalization.
"""
function clamp_site_fractions!(
    y::AbstractMatrix;
    min_val::Float64 = 1e-10,
    normalize::Bool = true,
)
    n_sub = size(y, 1)

    for s in 1:n_sub
        # Iterate to ensure min_val holds after normalization
        for _ in 1:10
            needs_clamp = false

            # Clamp to minimum
            for j in axes(y, 2)
                if y[s, j] < min_val
                    y[s, j] = min_val
                    needs_clamp = true
                end
            end

            # Normalize
            if normalize
                total = sum(y[s, :])
                if total > 0
                    y[s, :] ./= total
                end
            end

            # If no clamping was needed, constraint is satisfied
            if !needs_clamp
                break
            end
        end
    end

    return y
end

# Import LinearAlgebra for dot product
using LinearAlgebra: dot
