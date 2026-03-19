# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Gibbs energy calculation module

# =============================================================================
# Inden-Hillert-Jarl magnetic model
# =============================================================================

"""
    inden_hillert_g(tau, p)

Compute the g(τ) function in the Inden-Hillert-Jarl model.

  - τ < 1: Ferromagnetic region (below Curie temperature)
  - τ ≥ 1: Paramagnetic region (above Curie temperature)

# Arguments

  - `tau`: Reduced temperature T/Tc
  - `p`: Structure factor (0.40 for BCC, 0.28 for FCC/HCP)

# References

  - Hillert & Jarl (1978), Calphad 2(3), 227-238
  - Inden (1976), CALPHAD V, Chapter IX
"""
function inden_hillert_g(tau, p)
    A = 518 / 1125 + (11692 / 15975) * (1 / p - 1)
    if tau < 1
        term1 = 79 / (140 * p)
        term2 = (474 / 497) * (1 / p - 1) *
                (tau^3 / 6 + tau^9 / 135 + tau^15 / 600)
        g = 1 - (term1 * tau^(-1) + term2) / A
    else
        g = -(tau^(-5) / 10 + tau^(-15) / 315 +
              tau^(-25) / 1500) / A
    end
    return g
end

"""
    evaluate_magnetic_property(phase, T, y, db, symbol) -> Real

Evaluate a magnetic property (TC or BMAGN) as a composition-weighted
value using the same Redlich-Kister mixing as Gibbs energy.

# Arguments

  - `symbol`: `'T'` for Curie temperature, `'B'` for Bohr magneton
"""
function evaluate_magnetic_property(
    phase::Phase, T::TT, y::AbstractMatrix{Y},
    db::Database, P::Real, symbol::Char,
) where {TT <: Real, Y <: Real}
    RT = promote_type(TT, Y, Float64)
    val = zero(RT)

    for param in phase.parameters
        param.symbol == symbol || continue

        # Determine if endmember or interaction
        n_multi = count(c -> length(c) >= 2, param.constituents)

        if n_multi == 0
            # Endmember: product of site fractions * value
            y_prod = one(RT)
            for (s, constituents) in enumerate(param.constituents)
                if s > length(phase.constituents)
                    break
                end
                for spec in constituents
                    idx = findfirst(
                        ==(spec), phase.constituents[s],
                    )
                    if !isnothing(idx)
                        y_prod *= y[s, idx]
                    end
                end
            end
            if y_prod > 1e-30
                pval = evaluate_parameter(param, T, db, P)
                val += y_prod * pval
            end
        else
            # Interaction: Redlich-Kister
            contribution = one(RT)
            interaction_sublattice = 0
            interaction_species = String[]

            for (s, constituents) in enumerate(param.constituents)
                if s > length(phase.constituents)
                    break
                end
                if length(constituents) == 1
                    spec = constituents[1]
                    idx = findfirst(
                        ==(spec), phase.constituents[s],
                    )
                    if !isnothing(idx)
                        contribution *= y[s, idx]
                    else
                        contribution = zero(RT)
                        break
                    end
                elseif length(constituents) == 2
                    interaction_sublattice = s
                    interaction_species = constituents
                end
            end

            if contribution < 1e-30 ||
               interaction_sublattice == 0
                continue
            end

            spec1, spec2 = interaction_species
            idx1 = findfirst(
                ==(spec1),
                phase.constituents[interaction_sublattice],
            )
            idx2 = findfirst(
                ==(spec2),
                phase.constituents[interaction_sublattice],
            )
            if isnothing(idx1) || isnothing(idx2)
                continue
            end

            y1 = y[interaction_sublattice, idx1]
            y2 = y[interaction_sublattice, idx2]
            if y1 < 1e-30 || y2 < 1e-30
                continue
            end

            pval = evaluate_parameter(param, T, db, P)
            rk = y1 * y2 * (y1 - y2)^param.order * pval
            val += contribution * rk
        end
    end

    return val
end

"""
    calculate_magnetic_energy(phase, T, y, db; P=1e5) -> Real

Calculate the magnetic contribution to Gibbs energy using the
Inden-Hillert-Jarl model.

Returns 0 for phases without a magnetic model.

# Formula

G_mag = R * T * ln(β + 1) * g(T / Tc)

where Tc and β are composition-weighted averages of TC and BMAGN
parameters, and g(τ) is the Inden-Hillert-Jarl function.
"""
function calculate_magnetic_energy(
    phase::Phase, T::TT, y::AbstractMatrix{Y},
    db::Database; P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    RT = promote_type(TT, Y, Float64)

    mag = get(db.magnetic_models, uppercase(phase.name), nothing)
    isnothing(mag) && return zero(RT)

    Tc = evaluate_magnetic_property(phase, T, y, db, P, 'T')
    beta = evaluate_magnetic_property(
        phase, T, y, db, P, 'B',
    )

    # Anti-ferromagnetic handling
    if Tc < 0
        Tc = Tc / mag.afm_factor
        beta = beta / mag.afm_factor
    end

    # Guard: no magnetic contribution if Tc or beta is negligible
    (abs(Tc) < 1e-10 || abs(beta) < 1e-10) && return zero(RT)

    tau = T / Tc
    g = inden_hillert_g(tau, mag.p)
    return R * T * log(abs(beta) + 1) * g
end

"""
    calculate_gibbs_energy(phase::Phase, T::Real, y::AbstractMatrix, db::Database;
                          P::Real=1e5) -> Float64

Calculate the molar Gibbs energy of a phase at temperature T and pressure P
with site fractions y[sublattice, constituent].

The total Gibbs energy is:
G = G_ref + G_ideal + G_excess + G_magnetic

where:

  - G_ref: Reference state contribution (endmember energies)
  - G_ideal: Ideal mixing entropy contribution
  - G_excess: Excess Gibbs energy (Redlich-Kister)
  - G_magnetic: Magnetic contribution (if applicable)

# Examples

```julia
db = read_tdb(joinpath(pkgdir(OpenCALPHAD), "test", "data", "agcu.TDB"))
phase = get_phase(db, "FCC_A1")
y = make_y_matrix(phase, 0.3)
G = calculate_gibbs_energy(phase, 1000.0, y, db)
```
"""
function calculate_gibbs_energy(
    phase::Phase, T::TT, y::AbstractMatrix{Y}, db::Database;
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    RT = promote_type(TT, Y, Float64)
    G_ref = calculate_reference_energy(phase, T, y, db, P)::RT
    G_ideal = calculate_ideal_mixing(phase, T, y)::RT
    G_excess = calculate_excess_energy(phase, T, y, db, P)::RT
    G_mag = calculate_magnetic_energy(phase, T, y, db; P = P)::RT

    return G_ref + G_ideal + G_excess + G_mag
end

"""
    calculate_reference_energy(phase::Phase, T::Real, y::AbstractMatrix,
                              db::Database, P::Real) -> Float64

Calculate the reference state (endmember) contribution to Gibbs energy.
Uses the Compound Energy Formalism (CEF).

G_ref = Σ (∏ᵢ yᵢⱼ) * G^endmember

where the sum runs over all endmembers (combinations of constituents).
"""
function calculate_reference_energy(
    phase::Phase, T::TT, y::AbstractMatrix{Y},
    db::Database, P::Real,
) where {TT <: Real, Y <: Real}
    RT = promote_type(TT, Y, Float64)
    G_ref = zero(RT)

    # Find all G parameters that are endmembers (single species per sublattice)
    # Exclude interaction parameters (multiple species in any sublattice)
    for param in phase.parameters
        # Skip non-endmember parameters
        param.symbol == 'G' || continue
        all(length(c) == 1 for c in param.constituents) || continue

        # Calculate product of site fractions for this endmember
        y_prod = one(RT)
        for (s, constituents) in enumerate(param.constituents)
            if s > length(phase.constituents)
                break
            end
            for spec in constituents
                # Find index of this species in the phase constitution
                idx = findfirst(==(spec), phase.constituents[s])
                if !isnothing(idx)
                    y_prod *= y[s, idx]
                end
            end
        end

        # Evaluate the parameter and add contribution
        if y_prod > 1e-30
            G_param = evaluate_parameter(param, T, db, P)
            G_ref += y_prod * G_param
        end
    end

    return G_ref
end

"""
    calculate_ideal_mixing(phase::Phase, T::Real, y::AbstractMatrix) -> Float64

Calculate the ideal mixing contribution: R*T*Σᵢ(aᵢ*Σⱼ(yᵢⱼ*ln(yᵢⱼ)))
where aᵢ is the site multiplicity of sublattice i.
"""
function calculate_ideal_mixing(
    phase::Phase, T::TT, y::AbstractMatrix{Y},
) where {TT <: Real, Y <: Real}
    RT = promote_type(TT, Y, Float64)
    G_ideal = zero(RT)

    for (s, a) in enumerate(phase.sites)
        for j in axes(y, 2)
            yval = y[s, j]
            if yval > 1e-30  # Avoid log(0)
                G_ideal += a * yval * log(yval)
            end
        end
    end

    return R * T * G_ideal
end

"""
    calculate_excess_energy(phase::Phase, T::Real, y::AbstractMatrix,
                           db::Database, P::Real) -> Float64

Calculate the excess Gibbs energy using Redlich-Kister polynomial.
For binary interaction on sublattice s:
G_excess = Σₛ yᵢ⁽ˢ⁾ yⱼ⁽ˢ⁾ Σᵥ Lᵥ (yᵢ⁽ˢ⁾ - yⱼ⁽ˢ⁾)^v

In CEF/TDB format, interaction parameters are written as either:

  - PARAMETER L(PHASE,A,B:C;order) - explicit L parameters
  - PARAMETER G(PHASE,A,B:C;order) - G parameters with multiple species (interaction)
"""
function calculate_excess_energy(
    phase::Phase, T::TT, y::AbstractMatrix{Y},
    db::Database, P::Real,
) where {TT <: Real, Y <: Real}
    RT = promote_type(TT, Y, Float64)
    G_excess = zero(RT)

    for param in phase.parameters
        # Check if this is an interaction parameter
        is_interaction = false
        if param.symbol == 'L'
            is_interaction = true
        elseif param.symbol == 'G'
            for constituents in param.constituents
                if length(constituents) >= 2
                    is_interaction = true
                    break
                end
            end
        end
        is_interaction || continue

        # Calculate contribution from this interaction parameter
        # First, find which sublattice has the interaction (multiple species)
        contribution = one(RT)
        interaction_sublattice = 0
        interaction_species = String[]

        for (s, constituents) in enumerate(param.constituents)
            if s > length(phase.constituents)
                break
            end

            if length(constituents) == 1
                # Single species - just multiply by site fraction
                spec = constituents[1]
                idx = findfirst(==(spec), phase.constituents[s])
                if !isnothing(idx)
                    contribution *= y[s, idx]
                else
                    contribution = zero(RT)
                    break
                end
            elseif length(constituents) == 2
                # Binary interaction
                interaction_sublattice = s
                interaction_species = constituents
            end
        end

        if contribution < 1e-30 || interaction_sublattice == 0
            continue
        end

        # Get site fractions for interacting species
        spec1, spec2 = interaction_species
        idx1 = findfirst(==(spec1), phase.constituents[interaction_sublattice])
        idx2 = findfirst(==(spec2), phase.constituents[interaction_sublattice])

        if isnothing(idx1) || isnothing(idx2)
            continue
        end

        y1 = y[interaction_sublattice, idx1]
        y2 = y[interaction_sublattice, idx2]

        if y1 < 1e-30 || y2 < 1e-30
            continue
        end

        # Redlich-Kister: y1 * y2 * (y1 - y2)^order * L
        L = evaluate_parameter(param, T, db, P)
        rk_term = y1 * y2 * (y1 - y2)^param.order * L
        G_excess += contribution * rk_term
    end

    return G_excess
end

"""
    gibbs_energy_with_derivatives(phase::Phase, T::Real, y::AbstractMatrix,
                                  db::Database; P::Real=1e5)

Calculate Gibbs energy and its derivatives using automatic differentiation.
Returns (G, dG/dy, d²G/dy²).
"""
function gibbs_energy_with_derivatives(phase::Phase, T::Real, y::AbstractMatrix,
    db::Database; P::Real = 1e5)
    # Use ForwardDiff for automatic differentiation
    y_flat = vec(y)

    function G_func(y_vec)
        y_mat = reshape(y_vec, size(y))
        return calculate_gibbs_energy(phase, T, y_mat, db; P = P)
    end

    G = G_func(y_flat)
    dG = ForwardDiff.gradient(G_func, y_flat)
    d2G = ForwardDiff.hessian(G_func, y_flat)

    return (G, reshape(dG, size(y)), d2G)
end

# =============================================================================
# Julia DSL: Gibbs energy calculation using JuliaParameter
# =============================================================================

"""
    calculate_gibbs_energy(phase::Phase, T::Real, y::AbstractMatrix; P::Real=1e5) -> Float64

Calculate the molar Gibbs energy using Julia function parameters (DSL mode).
This method is used when the phase has `julia_parameters` defined.

No Database argument is needed since all functions are embedded in the phase.
"""
function calculate_gibbs_energy(
    phase::Phase, T::TT, y::AbstractMatrix{Y}; P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    if isempty(phase.julia_parameters)
        error(
            "Phase $(phase.name) has no Julia parameters. Use calculate_gibbs_energy(phase, T, y, db) with a Database.",
        )
    end

    RT = promote_type(TT, Y, Float64)
    G_ref = calculate_reference_energy_julia(phase, T, y)::RT
    G_ideal = calculate_ideal_mixing(phase, T, y)::RT
    G_excess = calculate_excess_energy_julia(phase, T, y)::RT

    return G_ref + G_ideal + G_excess
end

"""
    calculate_reference_energy_julia(phase::Phase, T::Real, y::AbstractMatrix) -> Float64

Calculate the reference state contribution using Julia function parameters.
"""
function calculate_reference_energy_julia(
    phase::Phase, T::TT, y::AbstractMatrix{Y},
) where {TT <: Real, Y <: Real}
    RT = promote_type(TT, Y, Float64)
    G_ref = zero(RT)

    for param in phase.julia_parameters
        # Skip non-endmember parameters
        param.symbol == 'G' || continue
        is_endmember(param) || continue

        # Calculate product of site fractions for this endmember
        y_prod = one(RT)
        for (s, constituents) in enumerate(param.constituents)
            if s > length(phase.constituents)
                break
            end
            for spec in constituents
                spec_str = String(spec)
                idx = findfirst(==(spec_str), phase.constituents[s])
                if !isnothing(idx)
                    y_prod *= y[s, idx]
                end
            end
        end

        # Evaluate the Julia function and add contribution
        # Use function barrier to isolate type-unstable call
        if y_prod > 1e-30
            G_param = _eval_julia_param(param.func, T)
            G_ref += y_prod * G_param
        end
    end

    return G_ref
end

#=
    _eval_julia_param(f::Function, T) -> Number

Evaluate a Julia DSL parameter function `f` at temperature `T`.
Marked @noinline as a function barrier to isolate the type-unstable call
and support ForwardDiff Dual numbers.
=#
@noinline _eval_julia_param(f::Function, T) = f(T)

"""
    calculate_excess_energy_julia(phase::Phase, T::Real, y::AbstractMatrix) -> Float64

Calculate the excess Gibbs energy using Julia function parameters.
Supports binary and higher-order interactions via Redlich-Kister or Muggianu formalism.
"""
function calculate_excess_energy_julia(
    phase::Phase, T::TT, y::AbstractMatrix{Y},
) where {TT <: Real, Y <: Real}
    RT = promote_type(TT, Y, Float64)
    G_excess = zero(RT)

    for param in phase.julia_parameters
        # Skip non-interaction parameters
        param.symbol == 'L' || continue

        # Find which sublattice has the interaction
        contribution = one(RT)
        interaction_sublattice = 0
        interaction_species = Symbol[]

        for (s, constituents) in enumerate(param.constituents)
            if s > length(phase.constituents)
                break
            end

            if length(constituents) == 1
                # Single species - multiply by site fraction
                spec = constituents[1]
                spec_str = String(spec)
                idx = findfirst(==(spec_str), phase.constituents[s])
                if !isnothing(idx)
                    contribution *= y[s, idx]
                else
                    contribution = zero(RT)
                    break
                end
            elseif length(constituents) == 2
                # Binary interaction
                interaction_sublattice = s
                interaction_species = constituents
            elseif length(constituents) >= 3
                # Ternary or higher interaction
                interaction_sublattice = s
                interaction_species = constituents
            end
        end

        if contribution < 1e-30 || interaction_sublattice == 0
            continue
        end

        # Calculate interaction term based on number of species
        if length(interaction_species) == 2
            # Binary Redlich-Kister: y1 * y2 * (y1 - y2)^order * L
            spec1, spec2 = interaction_species
            idx1 = findfirst(==(String(spec1)), phase.constituents[interaction_sublattice])
            idx2 = findfirst(==(String(spec2)), phase.constituents[interaction_sublattice])

            if isnothing(idx1) || isnothing(idx2)
                continue
            end

            y1 = y[interaction_sublattice, idx1]
            y2 = y[interaction_sublattice, idx2]

            if y1 < 1e-30 || y2 < 1e-30
                continue
            end

            L = _eval_julia_param(param.func, T)
            rk_term = y1 * y2 * (y1 - y2)^param.order * L
            G_excess += contribution * rk_term

        elseif length(interaction_species) >= 3
            # Ternary interaction: y1 * y2 * y3 * L (Muggianu)
            # Only order=0 for ternary
            if param.order != 0
                continue
            end

            y_vals = RT[]
            valid = true
            for spec in interaction_species
                idx =
                    findfirst(==(String(spec)), phase.constituents[interaction_sublattice])
                if isnothing(idx) || y[interaction_sublattice, idx] < 1e-30
                    valid = false
                    break
                end
                push!(y_vals, y[interaction_sublattice, idx])
            end

            if !valid
                continue
            end

            L = _eval_julia_param(param.func, T)
            # Product of all site fractions
            y_product = prod(y_vals)
            G_excess += contribution * y_product * L
        end
    end

    return G_excess
end

# =============================================================================
# Fallback methods for type errors
# =============================================================================

"""
    calculate_gibbs_energy(phase, T, y, db; P=1e5)

Fallback method for invalid phase type.
"""
function calculate_gibbs_energy(phase, T, y, db; kwargs...)
    error("Invalid phase type: $(typeof(phase)). Expected Phase.")
end

"""
    calculate_gibbs_energy(phase, T, y; P=1e5)

Fallback method for invalid phase type (DSL version).
"""
function calculate_gibbs_energy(phase, T, y; kwargs...)
    error("Invalid phase type: $(typeof(phase)). Expected Phase.")
end
