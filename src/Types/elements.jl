# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Element type definitions

"""
    Element

Represents a chemical element in the thermodynamic database.

# Fields

  - `name::String`: Element symbol (e.g., "AG", "CU")
  - `reference_state::String`: Reference phase (e.g., "FCC_A1")
  - `mass::Float64`: Atomic mass [g/mol]
  - `H298::Float64`: Enthalpy at 298.15 K [J/mol]
  - `S298::Float64`: Entropy at 298.15 K [J/(mol·K)]
"""
struct Element
    name::String
    reference_state::String
    mass::Float64
    H298::Float64
    S298::Float64
end

function Element(name::AbstractString, reference_state::AbstractString, mass::Real)
    return Element(String(name), String(reference_state), Float64(mass), 0.0, 0.0)
end

function Element(name::AbstractString, reference_state::AbstractString,
    mass::Real, H298::Real, S298::Real)
    return Element(String(name), String(reference_state),
        Float64(mass), Float64(H298), Float64(S298))
end

function Base.show(io::IO, elem::Element)
    print(io, "Element($(elem.name), ref=$(elem.reference_state), mass=$(elem.mass))")
end

# Special electron element (for charged species)
const ELECTRON = Element("/-", "ELECTRON_GAS", 0.0, 0.0, 0.0)

# Vacancy pseudo-element
const VACANCY = Element("VA", "VACUUM", 0.0, 0.0, 0.0)
