# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Phase Field Integration Utilities

"""
    driving_force(db::Database, T::Real, x::Real, solid_phase::String, liquid_phase::String;
                  P::Real=1e5) -> Float64

Calculate the thermodynamic driving force for solidification.

ΔG = G_solid(T, x) - G_liquid(T, x)

# Arguments

  - `db`: Database containing phase definitions
  - `T`: Temperature [K]
  - `x`: Mole fraction of second component (0 to 1)
  - `solid_phase`: Name of solid phase (e.g., "FCC_A1")
  - `liquid_phase`: Name of liquid phase (e.g., "LIQUID")
  - `P`: Pressure [Pa]

# Returns

  - ΔG [J/mol]: Negative means solidification is favorable

# Example

```julia
db = read_tdb("agcu.TDB")
ΔG = driving_force(db, 1000.0, 0.3, "FCC_A1", "LIQUID")
```
"""
function driving_force(
    db::Database,
    T::Real,
    x::Real,
    solid_phase::String,
    liquid_phase::String;
    P::Real = 1e5,
)
    solid = get_phase(db, solid_phase)
    liquid = get_phase(db, liquid_phase)

    y = _composition_to_site_fractions(solid, x)
    G_s = calculate_gibbs_energy(solid, T, y, db; P = P)

    y_l = _composition_to_site_fractions(liquid, x)
    G_l = calculate_gibbs_energy(liquid, T, y_l, db; P = P)

    return G_s - G_l
end

"""
    chemical_potential(phase::Phase, T::Real, x::Real, db::Database;
                       P::Real=1e5) -> Tuple{Float64, Float64}

Calculate chemical potentials of both components.

μ_i = ∂G/∂n_i at constant T, P, n_j

# Arguments

  - `phase`: Phase object
  - `T`: Temperature [K]
  - `x`: Mole fraction of second component
  - `db`: Database
  - `P`: Pressure [Pa]

# Returns

  - (μ₁, μ₂): Chemical potentials [J/mol] of components 1 and 2

# Example

```julia
fcc = get_phase(db, "FCC_A1")
μ_Ag, μ_Cu = chemical_potential(fcc, 1000.0, 0.3, db)
```
"""
function chemical_potential(
    phase::Phase,
    T::Real,
    x::Real,
    db::Database;
    P::Real = 1e5,
)
    y = _composition_to_site_fractions(phase, x)
    G, dG, _ = gibbs_energy_with_derivatives(phase, T, y, db; P = P)

    # For binary: μ₁ = G - x·(∂G/∂x), μ₂ = G + (1-x)·(∂G/∂x)
    # dG is with respect to site fractions, need to convert
    dGdx = dG[1, 2] - dG[1, 1]  # ∂G/∂x = ∂G/∂y₂ - ∂G/∂y₁

    μ₁ = G - x * dGdx
    μ₂ = G + (1 - x) * dGdx

    return (μ₁, μ₂)
end

"""
    chemical_potential(phase::Phase, T::Real, x::Real; P::Real=1e5) -> Tuple{Float64, Float64}

Calculate chemical potentials using Julia DSL parameters (no Database required).
"""
function chemical_potential(
    phase::Phase,
    T::Real,
    x::Real;
    P::Real = 1e5,
)
    if !has_julia_parameters(phase)
        error(
            "Phase $(phase.name) has no Julia parameters. Use chemical_potential(phase, T, x, db).",
        )
    end

    y = _composition_to_site_fractions(phase, x)

    # Use DifferentiationInterface to get derivatives
    function G_func(x_val)
        y_local = _composition_to_site_fractions(phase, x_val)
        return calculate_gibbs_energy(phase, T, y_local; P = P)
    end

    G = G_func(x)
    dGdx = DI.derivative(G_func, DI.AutoForwardDiff(), x)

    μ₁ = G - x * dGdx
    μ₂ = G + (1 - x) * dGdx

    return (μ₁, μ₂)
end

"""
    diffusion_potential(phase::Phase, T::Real, x::Real, db::Database;
                        P::Real=1e5) -> Float64

Calculate the diffusion potential (second derivative of G).

∂²G/∂x² - used in Cahn-Hilliard equation for diffusion coefficient.

# Returns

  - d²G/dx² [J/mol]
"""
function diffusion_potential(
    phase::Phase,
    T::Real,
    x::Real,
    db::Database;
    P::Real = 1e5,
)
    function G_func(x_val)
        y = _composition_to_site_fractions(phase, x_val)
        return calculate_gibbs_energy(phase, T, y, db; P = P)
    end

    # Second derivative using DifferentiationInterface
    d2Gdx2 = DI.second_derivative(G_func, DI.AutoForwardDiff(), x)

    return d2Gdx2
end

"""
    phase_field_params(db::Database, T::Real, x::Real,
                       solid_phase::String, liquid_phase::String;
                       P::Real=1e5) -> NamedTuple

Calculate all parameters needed for phase field simulation.

# Returns

NamedTuple with:

  - `ΔG`: Driving force [J/mol]
  - `μ_solid`: Chemical potentials in solid (μ₁, μ₂)
  - `μ_liquid`: Chemical potentials in liquid (μ₁, μ₂)
  - `d2G_solid`: Diffusion potential in solid
  - `d2G_liquid`: Diffusion potential in liquid
"""
function phase_field_params(
    db::Database,
    T::Real,
    x::Real,
    solid_phase::String,
    liquid_phase::String;
    P::Real = 1e5,
)
    solid = get_phase(db, solid_phase)
    liquid = get_phase(db, liquid_phase)

    ΔG = driving_force(db, T, x, solid_phase, liquid_phase; P = P)
    μ_s = chemical_potential(solid, T, x, db; P = P)
    μ_l = chemical_potential(liquid, T, x, db; P = P)
    d2G_s = diffusion_potential(solid, T, x, db; P = P)
    d2G_l = diffusion_potential(liquid, T, x, db; P = P)

    return (
        ΔG = ΔG,
        μ_solid = μ_s,
        μ_liquid = μ_l,
        d2G_solid = d2G_s,
        d2G_liquid = d2G_l,
    )
end

# =============================================================================
# Helper functions
# =============================================================================

"""
Convert mole fraction x to site fraction matrix for binary phase.
"""
function _composition_to_site_fractions(phase::Phase, x::Real)
    n_sub = length(phase.sites)
    n_const = maximum(length.(phase.constituents))

    y = zeros(typeof(x), n_sub, n_const)
    y[1, 1] = 1 - x  # First component
    y[1, 2] = x      # Second component

    # Fill other sublattices with VA
    for s in 2:n_sub
        va_idx = find_va_index(phase.constituents[s])
        y[s, va_idx] = one(x)
    end

    return y
end
