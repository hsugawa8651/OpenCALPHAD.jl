# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Equilibrium result type definitions

"""
    PhaseResult

Results for a single phase in equilibrium.

# Fields

  - `name::String`: Phase name
  - `status::Symbol`: Phase status (:stable, :metastable, :suspended)
  - `amount::Float64`: Amount of phase [mol]
  - `site_fractions::Matrix{Float64}`: Site fractions y[sublattice, constituent]
  - `mole_fractions::Vector{Float64}`: Mole fractions of components in phase
  - `gibbs_energy::Float64`: Molar Gibbs energy [J/mol]
  - `driving_force::Float64`: Driving force for precipitation [J/mol]
"""
mutable struct PhaseResult
    name::String
    status::Symbol
    amount::Float64
    site_fractions::Matrix{Float64}
    mole_fractions::Vector{Float64}
    gibbs_energy::Float64
    driving_force::Float64
end

function PhaseResult(name::AbstractString)
    return PhaseResult(
        String(name),
        :suspended,
        0.0,
        Matrix{Float64}(undef, 0, 0),
        Float64[],
        0.0,
        0.0,
    )
end

function Base.show(io::IO, pr::PhaseResult)
    status_str = pr.status == :stable ? "E" : (pr.status == :metastable ? "M" : "S")
    print(io, "$(pr.name) [$status_str]: $(pr.amount) mol")
end

"""
    ComponentResult

Results for a single component in equilibrium.

# Fields

  - `name::String`: Component name
  - `moles::Float64`: Total moles
  - `mole_fraction::Float64`: Overall mole fraction
  - `chemical_potential::Float64`: Chemical potential [J/mol]
  - `activity::Float64`: Activity (relative to reference state)
"""
struct ComponentResult
    name::String
    moles::Float64
    mole_fraction::Float64
    chemical_potential::Float64
    activity::Float64
end

function Base.show(io::IO, cr::ComponentResult)
    @printf(io, "%s: n=%.4e, x=%.4f, μ/RT=%.4f, a=%.4e",
        cr.name, cr.moles, cr.mole_fraction,
        cr.chemical_potential / (R * 1000), cr.activity)
end

"""
    Equilibrium

Complete equilibrium calculation result.

# Fields

  - `conditions::ConditionSet`: Input conditions
  - `components::Vector{ComponentResult}`: Component results
  - `phases::Vector{PhaseResult}`: Phase results
  - `gibbs_energy::Float64`: Total Gibbs energy [J]
  - `enthalpy::Float64`: Total enthalpy [J]
  - `entropy::Float64`: Total entropy [J/K]
  - `iterations::Int`: Number of iterations to converge
  - `converged::Bool`: Whether calculation converged
"""
mutable struct Equilibrium
    conditions::ConditionSet
    components::Vector{ComponentResult}
    phases::Vector{PhaseResult}
    gibbs_energy::Float64
    enthalpy::Float64
    entropy::Float64
    iterations::Int
    converged::Bool
end

function Equilibrium(conditions::ConditionSet)
    return Equilibrium(
        conditions,
        ComponentResult[],
        PhaseResult[],
        0.0,
        0.0,
        0.0,
        0,
        false,
    )
end

function Base.show(io::IO, eq::Equilibrium)
    println(io, "Equilibrium Result:")
    println(io, "  $(eq.conditions)")
    println(io, "  Converged: $(eq.converged) ($(eq.iterations) iterations)")
    println(io, "  G = $(eq.gibbs_energy) J")
    println(io, "  Components:")
    for comp in eq.components
        println(io, "    $comp")
    end
    println(io, "  Phases:")
    for phase in eq.phases
        if phase.status == :stable
            println(io, "    $phase")
        end
    end
end

"""
    stable_phases(eq::Equilibrium) -> Vector{PhaseResult}

Return only the stable phases from equilibrium result.
"""
function stable_phases(eq::Equilibrium)
    return filter(p -> p.status == :stable, eq.phases)
end

"""
    n_stable_phases(eq::Equilibrium) -> Int

Return the number of stable phases.
"""
function n_stable_phases(eq::Equilibrium)
    return count(p -> p.status == :stable, eq.phases)
end
