# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Thermodynamic property functions (S, H, mixing properties)

"""
    calculate_entropy(phase::Phase, T::Real, y::AbstractMatrix, db::Database;
                      P::Real=1e5) -> Float64

Calculate the molar entropy S = -∂G/∂T using automatic differentiation.

Currently only supports binary (2-component) systems.

# Arguments

  - `phase`: Phase structure from database
  - `T`: Temperature [K]
  - `y`: Site fraction matrix y[sublattice, constituent]
  - `db`: Thermodynamic database
  - `P`: Pressure [Pa] (default: 1e5)

# Returns

  - Molar entropy [J/(mol·K)]
"""
function calculate_entropy(
    phase::Phase, T::TT, y::AbstractMatrix{Y}, db::Database;
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    # Check component count (currently only binary supported)
    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_entropy: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    # S = -∂G/∂T using ForwardDiff
    dG_dT = ForwardDiff.derivative(
        t -> calculate_gibbs_energy(phase, t, y, db; P = P),
        T,
    )
    return -dG_dT
end

"""
    calculate_enthalpy(phase::Phase, T::Real, y::AbstractMatrix, db::Database;
                       P::Real=1e5) -> Float64

Calculate the molar enthalpy H = G + T*S.

Currently only supports binary (2-component) systems.

# Arguments

  - `phase`: Phase structure from database
  - `T`: Temperature [K]
  - `y`: Site fraction matrix y[sublattice, constituent]
  - `db`: Thermodynamic database
  - `P`: Pressure [Pa] (default: 1e5)

# Returns

  - Molar enthalpy [J/mol]
"""
function calculate_enthalpy(
    phase::Phase, T::TT, y::AbstractMatrix{Y}, db::Database;
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    # Check component count
    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_enthalpy: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    G = calculate_gibbs_energy(phase, T, y, db; P = P)
    S = calculate_entropy(phase, T, y, db; P = P)
    return G + T * S
end

"""
    calculate_mixing_gibbs(phase::Phase, T::Real, y::AbstractMatrix, db::Database;
                           P::Real=1e5) -> Float64

Calculate the mixing Gibbs energy: G_mix = G - Σ x_i * G_i^pure

Currently only supports binary (2-component) systems.

# Arguments

  - `phase`: Phase structure from database
  - `T`: Temperature [K]
  - `y`: Site fraction matrix y[sublattice, constituent]
  - `db`: Thermodynamic database
  - `P`: Pressure [Pa] (default: 1e5)

# Returns

  - Mixing Gibbs energy [J/mol]
"""
function calculate_mixing_gibbs(
    phase::Phase, T::TT, y::AbstractMatrix{Y}, db::Database;
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    # Check component count
    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_mixing_gibbs: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    # Total Gibbs energy
    G = calculate_gibbs_energy(phase, T, y, db; P = P)

    # Mole fractions from first sublattice (binary system)
    x1 = y[1, 1]
    x2 = y[1, 2]

    # Pure element Gibbs energies
    G1_pure = _calculate_pure_gibbs(phase, T, 1, db; P = P)
    G2_pure = _calculate_pure_gibbs(phase, T, 2, db; P = P)

    # Mixing Gibbs energy
    G_ref = x1 * G1_pure + x2 * G2_pure
    return G - G_ref
end

"""
    calculate_mixing_enthalpy(phase::Phase, T::Real, y::AbstractMatrix, db::Database;
                              P::Real=1e5) -> Float64

Calculate the mixing enthalpy: H_mix = H - Σ x_i * H_i^pure

Currently only supports binary (2-component) systems.

# Arguments

  - `phase`: Phase structure from database
  - `T`: Temperature [K]
  - `y`: Site fraction matrix y[sublattice, constituent]
  - `db`: Thermodynamic database
  - `P`: Pressure [Pa] (default: 1e5)

# Returns

  - Mixing enthalpy [J/mol]
"""
function calculate_mixing_enthalpy(
    phase::Phase, T::TT, y::AbstractMatrix{Y}, db::Database;
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    # Check component count
    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_mixing_enthalpy: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    # Total enthalpy
    H = calculate_enthalpy(phase, T, y, db; P = P)

    # Mole fractions from first sublattice (binary system)
    x1 = y[1, 1]
    x2 = y[1, 2]

    # Pure element enthalpies
    H1_pure = _calculate_pure_enthalpy(phase, T, 1, db; P = P)
    H2_pure = _calculate_pure_enthalpy(phase, T, 2, db; P = P)

    # Mixing enthalpy
    H_ref = x1 * H1_pure + x2 * H2_pure
    return H - H_ref
end

"""
    calculate_mixing_entropy(phase::Phase, T::Real, y::AbstractMatrix, db::Database;
                             P::Real=1e5) -> Float64

Calculate the mixing entropy: S_mix = S - Σ x_i * S_i^pure

Currently only supports binary (2-component) systems.

# Arguments

  - `phase`: Phase structure from database
  - `T`: Temperature [K]
  - `y`: Site fraction matrix y[sublattice, constituent]
  - `db`: Thermodynamic database
  - `P`: Pressure [Pa] (default: 1e5)

# Returns

  - Mixing entropy [J/(mol·K)]
"""
function calculate_mixing_entropy(
    phase::Phase, T::TT, y::AbstractMatrix{Y}, db::Database;
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    # Check component count
    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_mixing_entropy: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    # Total entropy
    S = calculate_entropy(phase, T, y, db; P = P)

    # Mole fractions from first sublattice (binary system)
    x1 = y[1, 1]
    x2 = y[1, 2]

    # Pure element entropies
    S1_pure = _calculate_pure_entropy(phase, T, 1, db; P = P)
    S2_pure = _calculate_pure_entropy(phase, T, 2, db; P = P)

    # Mixing entropy
    S_ref = x1 * S1_pure + x2 * S2_pure
    return S - S_ref
end

"""
    thermodynamic_properties(phase::Phase, T::Real, y::AbstractMatrix, db::Database;
                             P::Real=1e5) -> NamedTuple

Calculate all thermodynamic properties efficiently in one call.

Currently only supports binary (2-component) systems.

# Returns

NamedTuple with fields:

  - `G`: Gibbs energy [J/mol]
  - `H`: Enthalpy [J/mol]
  - `S`: Entropy [J/(mol·K)]
  - `G_mix`: Mixing Gibbs energy [J/mol]
  - `H_mix`: Mixing enthalpy [J/mol]
  - `S_mix`: Mixing entropy [J/(mol·K)]
"""
function thermodynamic_properties(
    phase::Phase, T::TT, y::AbstractMatrix{Y}, db::Database;
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    # Check component count
    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "thermodynamic_properties: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    # Calculate total properties (efficient: reuse G calculation for S)
    G = calculate_gibbs_energy(phase, T, y, db; P = P)
    dG_dT = ForwardDiff.derivative(
        t -> calculate_gibbs_energy(phase, t, y, db; P = P),
        T,
    )
    S = -dG_dT
    H = G + T * S

    # Mole fractions
    x1 = y[1, 1]
    x2 = y[1, 2]

    # Pure element properties (calculated once each)
    G1_pure = _calculate_pure_gibbs(phase, T, 1, db; P = P)
    G2_pure = _calculate_pure_gibbs(phase, T, 2, db; P = P)

    dG1_dT = ForwardDiff.derivative(
        t -> _calculate_pure_gibbs(phase, t, 1, db; P = P),
        T,
    )
    dG2_dT = ForwardDiff.derivative(
        t -> _calculate_pure_gibbs(phase, t, 2, db; P = P),
        T,
    )

    S1_pure = -dG1_dT
    S2_pure = -dG2_dT
    H1_pure = G1_pure + T * S1_pure
    H2_pure = G2_pure + T * S2_pure

    # Mixing properties
    G_ref = x1 * G1_pure + x2 * G2_pure
    H_ref = x1 * H1_pure + x2 * H2_pure
    S_ref = x1 * S1_pure + x2 * S2_pure

    G_mix = G - G_ref
    H_mix = H - H_ref
    S_mix = S - S_ref

    return (G = G, H = H, S = S, G_mix = G_mix, H_mix = H_mix, S_mix = S_mix)
end

# =============================================================================
# DSL versions (without Database argument)
# =============================================================================

"""
    calculate_entropy(phase::Phase, T::Real, y::AbstractMatrix; P::Real=1e5) -> Float64

DSL version: Calculate entropy using Julia-defined parameters (no database needed).
"""
function calculate_entropy(
    phase::Phase, T::TT, y::AbstractMatrix{Y};
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    if isempty(phase.julia_parameters)
        error(
            "Phase $(phase.name) has no Julia parameters. Use calculate_entropy(phase, T, y, db) with a Database.",
        )
    end

    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_entropy: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    dG_dT = ForwardDiff.derivative(
        t -> calculate_gibbs_energy(phase, t, y; P = P),
        T,
    )
    return -dG_dT
end

"""
    calculate_enthalpy(phase::Phase, T::Real, y::AbstractMatrix; P::Real=1e5) -> Float64

DSL version: Calculate enthalpy using Julia-defined parameters (no database needed).
"""
function calculate_enthalpy(
    phase::Phase, T::TT, y::AbstractMatrix{Y};
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    if isempty(phase.julia_parameters)
        error(
            "Phase $(phase.name) has no Julia parameters. Use calculate_enthalpy(phase, T, y, db) with a Database.",
        )
    end

    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_enthalpy: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    G = calculate_gibbs_energy(phase, T, y; P = P)
    S = calculate_entropy(phase, T, y; P = P)
    return G + T * S
end

"""
    calculate_mixing_gibbs(phase::Phase, T::Real, y::AbstractMatrix; P::Real=1e5) -> Float64

DSL version: Calculate mixing Gibbs energy using Julia-defined parameters.
"""
function calculate_mixing_gibbs(
    phase::Phase, T::TT, y::AbstractMatrix{Y};
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    if isempty(phase.julia_parameters)
        error(
            "Phase $(phase.name) has no Julia parameters. Use calculate_mixing_gibbs(phase, T, y, db) with a Database.",
        )
    end

    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_mixing_gibbs: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    G = calculate_gibbs_energy(phase, T, y; P = P)
    x1 = y[1, 1]
    x2 = y[1, 2]
    G1_pure = _calculate_pure_gibbs_dsl(phase, T, 1; P = P)
    G2_pure = _calculate_pure_gibbs_dsl(phase, T, 2; P = P)
    G_ref = x1 * G1_pure + x2 * G2_pure
    return G - G_ref
end

"""
    calculate_mixing_enthalpy(phase::Phase, T::Real, y::AbstractMatrix; P::Real=1e5) -> Float64

DSL version: Calculate mixing enthalpy using Julia-defined parameters.
"""
function calculate_mixing_enthalpy(
    phase::Phase, T::TT, y::AbstractMatrix{Y};
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    if isempty(phase.julia_parameters)
        error(
            "Phase $(phase.name) has no Julia parameters. Use calculate_mixing_enthalpy(phase, T, y, db) with a Database.",
        )
    end

    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_mixing_enthalpy: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    H = calculate_enthalpy(phase, T, y; P = P)
    x1 = y[1, 1]
    x2 = y[1, 2]
    H1_pure = _calculate_pure_enthalpy_dsl(phase, T, 1; P = P)
    H2_pure = _calculate_pure_enthalpy_dsl(phase, T, 2; P = P)
    H_ref = x1 * H1_pure + x2 * H2_pure
    return H - H_ref
end

"""
    calculate_mixing_entropy(phase::Phase, T::Real, y::AbstractMatrix; P::Real=1e5) -> Float64

DSL version: Calculate mixing entropy using Julia-defined parameters.
"""
function calculate_mixing_entropy(
    phase::Phase, T::TT, y::AbstractMatrix{Y};
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    if isempty(phase.julia_parameters)
        error(
            "Phase $(phase.name) has no Julia parameters. Use calculate_mixing_entropy(phase, T, y, db) with a Database.",
        )
    end

    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "calculate_mixing_entropy: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    S = calculate_entropy(phase, T, y; P = P)
    x1 = y[1, 1]
    x2 = y[1, 2]
    S1_pure = _calculate_pure_entropy_dsl(phase, T, 1; P = P)
    S2_pure = _calculate_pure_entropy_dsl(phase, T, 2; P = P)
    S_ref = x1 * S1_pure + x2 * S2_pure
    return S - S_ref
end

"""
    thermodynamic_properties(phase::Phase, T::Real, y::AbstractMatrix; P::Real=1e5) -> NamedTuple

DSL version: Calculate all thermodynamic properties using Julia-defined parameters.
"""
function thermodynamic_properties(
    phase::Phase, T::TT, y::AbstractMatrix{Y};
    P::Real = 1e5,
) where {TT <: Real, Y <: Real}
    if isempty(phase.julia_parameters)
        error(
            "Phase $(phase.name) has no Julia parameters. Use thermodynamic_properties(phase, T, y, db) with a Database.",
        )
    end

    n_components = length(phase.constituents[1])
    if n_components != 2
        throw(
            ArgumentError(
                "thermodynamic_properties: only binary (2-component) systems supported. " *
                "Got $n_components components.",
            ),
        )
    end

    G = calculate_gibbs_energy(phase, T, y; P = P)
    dG_dT = ForwardDiff.derivative(t -> calculate_gibbs_energy(phase, t, y; P = P), T)
    S = -dG_dT
    H = G + T * S

    x1 = y[1, 1]
    x2 = y[1, 2]

    G1_pure = _calculate_pure_gibbs_dsl(phase, T, 1; P = P)
    G2_pure = _calculate_pure_gibbs_dsl(phase, T, 2; P = P)
    dG1_dT = ForwardDiff.derivative(t -> _calculate_pure_gibbs_dsl(phase, t, 1; P = P), T)
    dG2_dT = ForwardDiff.derivative(t -> _calculate_pure_gibbs_dsl(phase, t, 2; P = P), T)
    S1_pure = -dG1_dT
    S2_pure = -dG2_dT
    H1_pure = G1_pure + T * S1_pure
    H2_pure = G2_pure + T * S2_pure

    G_ref = x1 * G1_pure + x2 * G2_pure
    H_ref = x1 * H1_pure + x2 * H2_pure
    S_ref = x1 * S1_pure + x2 * S2_pure

    G_mix = G - G_ref
    H_mix = H - H_ref
    S_mix = S - S_ref

    return (G = G, H = H, S = S, G_mix = G_mix, H_mix = H_mix, S_mix = S_mix)
end

# =============================================================================
# Internal helper functions
# =============================================================================

"""
    _make_pure_y(phase::Phase, component_idx::Int) -> Matrix{Float64}

Create site fraction matrix for pure component (all other components = 0).
"""
function _make_pure_y(phase::Phase, component_idx::Int)
    n_sub = length(phase.constituents)
    n_const = maximum(length(c) for c in phase.constituents)
    y = zeros(n_sub, n_const)

    # First sublattice: only specified component = 1
    y[1, component_idx] = 1.0

    # Other sublattices: VA = 1 (auto-detect VA position)
    for s in 2:n_sub
        va_idx = find_va_index(phase.constituents[s])
        y[s, va_idx] = 1.0
    end

    return y
end

"""
    _calculate_pure_gibbs(phase, T, component_idx, db; P=1e5) -> Float64

Calculate Gibbs energy of pure component in the phase.
"""
function _calculate_pure_gibbs(
    phase::Phase, T::Real, component_idx::Int, db::Database;
    P::Real = 1e5,
)
    y_pure = _make_pure_y(phase, component_idx)
    return calculate_gibbs_energy(phase, T, y_pure, db; P = P)
end

"""
    _calculate_pure_entropy(phase, T, component_idx, db; P=1e5) -> Float64

Calculate entropy of pure component in the phase.
"""
function _calculate_pure_entropy(
    phase::Phase, T::Real, component_idx::Int, db::Database;
    P::Real = 1e5,
)
    dG_dT = ForwardDiff.derivative(
        t -> _calculate_pure_gibbs(phase, t, component_idx, db; P = P),
        T,
    )
    return -dG_dT
end

"""
    _calculate_pure_enthalpy(phase, T, component_idx, db; P=1e5) -> Float64

Calculate enthalpy of pure component in the phase.
"""
function _calculate_pure_enthalpy(
    phase::Phase, T::Real, component_idx::Int, db::Database;
    P::Real = 1e5,
)
    G = _calculate_pure_gibbs(phase, T, component_idx, db; P = P)
    S = _calculate_pure_entropy(phase, T, component_idx, db; P = P)
    return G + T * S
end

# DSL versions of pure component helpers

"""
    _calculate_pure_gibbs_dsl(phase, T, component_idx; P=1e5) -> Float64

DSL version: Calculate Gibbs energy of pure component.
"""
function _calculate_pure_gibbs_dsl(phase::Phase, T::Real, component_idx::Int; P::Real = 1e5)
    y_pure = _make_pure_y(phase, component_idx)
    return calculate_gibbs_energy(phase, T, y_pure; P = P)
end

"""
    _calculate_pure_entropy_dsl(phase, T, component_idx; P=1e5) -> Float64

DSL version: Calculate entropy of pure component.
"""
function _calculate_pure_entropy_dsl(
    phase::Phase,
    T::Real,
    component_idx::Int;
    P::Real = 1e5,
)
    dG_dT = ForwardDiff.derivative(
        t -> _calculate_pure_gibbs_dsl(phase, t, component_idx; P = P),
        T,
    )
    return -dG_dT
end

"""
    _calculate_pure_enthalpy_dsl(phase, T, component_idx; P=1e5) -> Float64

DSL version: Calculate enthalpy of pure component.
"""
function _calculate_pure_enthalpy_dsl(
    phase::Phase,
    T::Real,
    component_idx::Int;
    P::Real = 1e5,
)
    G = _calculate_pure_gibbs_dsl(phase, T, component_idx; P = P)
    S = _calculate_pure_entropy_dsl(phase, T, component_idx; P = P)
    return G + T * S
end
