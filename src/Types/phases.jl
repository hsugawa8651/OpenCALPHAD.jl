# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Phase type definitions

"""
    Sublattice

Type alias for sublattice constituents (vector of species names).
"""
const Sublattice = Vector{String}

"""
    Constitution

Type alias for phase constitution (vector of sublattices).
"""
const Constitution = Vector{Sublattice}

"""
    Parameter

Represents a thermodynamic parameter (G, L, TC, BMAGN, etc.) in TDB format.

# Fields

  - `phase_name::String`: Phase this parameter belongs to
  - `symbol::Char`: Parameter type ('G', 'L', 'T', 'B', etc.)
  - `constituents::Constitution`: Species combination for this parameter
  - `order::Int`: Interaction order (for L parameters: 0, 1, 2, ...)
  - `temp_range::Tuple{Float64,Float64}`: Valid temperature range
  - `temp_breaks::Vector{Float64}`: Temperature breakpoints
  - `expressions::Vector{String}`: Expression strings
"""
struct Parameter
    phase_name::String
    symbol::Char
    constituents::Constitution
    order::Int
    temp_range::Tuple{Float64, Float64}
    temp_breaks::Vector{Float64}
    expressions::Vector{String}
end

function Parameter(phase_name::AbstractString,
    symbol::Char,
    constituents::Vector{<:Vector},
    order::Int,
    temp_breaks::Vector{<:Real},
    expressions::Vector{<:AbstractString})
    temp_range = (minimum(temp_breaks), maximum(temp_breaks))
    return Parameter(
        String(phase_name),
        symbol,
        [String.(c) for c in constituents],
        order,
        temp_range,
        Float64.(temp_breaks),
        String.(expressions),
    )
end

function Base.show(io::IO, p::Parameter)
    cons_str = join([join(s, ",") for s in p.constituents], ":")
    print(io, "Parameter($(p.symbol)($(p.phase_name),$(cons_str);$(p.order)))")
end

"""
    JuliaParameter

Represents a thermodynamic parameter defined as a Julia function.
This allows direct use of Julia code instead of TDB expression strings.

# Fields

  - `symbol::Char`: Parameter type ('G' for Gibbs, 'L' for interaction)
  - `constituents::Vector{Vector{Symbol}}`: Species per sublattice
  - `order::Int`: Interaction order (for L parameters: 0, 1, 2, ...)
  - `func::Function`: Julia function T -> Float64
"""
struct JuliaParameter
    symbol::Char
    constituents::Vector{Vector{Symbol}}
    order::Int
    func::Function
end

# Convenience constructors
function JuliaParameter(symbol::Char, constituents::Vector{Symbol}, func::Function)
    # Single sublattice: wrap in vector
    return JuliaParameter(symbol, [constituents], 0, func)
end

function JuliaParameter(
    symbol::Char,
    constituents::Vector{Symbol},
    order::Int,
    func::Function,
)
    return JuliaParameter(symbol, [constituents], order, func)
end

function Base.show(io::IO, p::JuliaParameter)
    cons_str = join([join(string.(s), ",") for s in p.constituents], ":")
    print(io, "JuliaParameter($(p.symbol)($(cons_str);$(p.order)))")
end

"""
    is_endmember(p::JuliaParameter) -> Bool

Check if the parameter is an endmember (single species per sublattice).
"""
is_endmember(p::JuliaParameter) = all(length(c) == 1 for c in p.constituents)

"""
    is_interaction(p::JuliaParameter) -> Bool

Check if the parameter is an interaction (multiple species in some sublattice).
"""
is_interaction(p::JuliaParameter) = any(length(c) >= 2 for c in p.constituents)

"""
    PhaseType

Enumeration of phase model types.
"""
@enum PhaseType begin
    PURE_SUBSTANCE      # Single component
    IDEAL_SOLUTION      # Ideal mixing
    REGULAR_SOLUTION    # Regular solution (Redlich-Kister)
    SUBLATTICE_MODEL    # Compound Energy Formalism (CEF)
    IONIC_LIQUID        # Two-sublattice ionic liquid
    AQUEOUS             # Aqueous solution
    GAS                 # Gas phase (ideal gas mixture)
end

"""
    Phase

Represents a thermodynamic phase.

# Fields

  - `name::String`: Phase identifier (e.g., "FCC_A1", "LIQUID")
  - `state::Char`: State indicator ('L', 'S', 'G')
  - `model_type::String`: Model type string from TDB
  - `sites::Vector{Float64}`: Site ratios for each sublattice
  - `constituents::Constitution`: Allowed species on each sublattice
  - `parameters::Vector{Parameter}`: Thermodynamic parameters (from TDB)
  - `julia_parameters::Vector{JuliaParameter}`: Parameters defined as Julia functions
"""
mutable struct Phase
    name::String
    state::Char
    model_type::String
    sites::Vector{Float64}
    constituents::Constitution
    parameters::Vector{Parameter}
    julia_parameters::Vector{JuliaParameter}
end

# Constructor with all fields
function Phase(name::AbstractString, state::Char, model_type::AbstractString,
    sites::Vector{<:Real}, constituents::Constitution,
    parameters::Vector{Parameter}, julia_parameters::Vector{JuliaParameter})
    return Phase(
        String(name),
        state,
        String(model_type),
        Float64.(sites),
        constituents,
        parameters,
        julia_parameters,
    )
end

# Backward compatible constructors (without julia_parameters)
function Phase(name::AbstractString, state::Char, sites::Vector{<:Real})
    return Phase(
        String(name),
        state,
        "%",  # Default model type
        Float64.(sites),
        Constitution(),
        Parameter[],
        JuliaParameter[],
    )
end

function Phase(name::AbstractString, state::Char, model_type::AbstractString,
    sites::Vector{<:Real})
    return Phase(
        String(name),
        state,
        String(model_type),
        Float64.(sites),
        Constitution(),
        Parameter[],
        JuliaParameter[],
    )
end

# DSL-style constructor: Phase(name, sites, constituents)
function Phase(name::AbstractString, sites::Vector{<:Real},
    constituents::Vector{Vector{Symbol}})
    const_str = [String.(c) for c in constituents]
    return Phase(
        String(name),
        'S',  # Default to solid
        "%",
        Float64.(sites),
        const_str,
        Parameter[],
        JuliaParameter[],
    )
end

function Base.show(io::IO, p::Phase)
    nsub = length(p.sites)
    npar = length(p.parameters) + length(p.julia_parameters)
    print(io, "Phase($(p.name), state=$(p.state), sublattices=$nsub, params=$npar)")
end

"""
    has_julia_parameters(phase::Phase) -> Bool

Check if the phase has any Julia function parameters defined.
"""
has_julia_parameters(phase::Phase) = !isempty(phase.julia_parameters)

"""
    num_sublattices(phase::Phase) -> Int

Return the number of sublattices in the phase.
"""
num_sublattices(phase::Phase) = length(phase.sites)

"""
    total_sites(phase::Phase) -> Float64

Return the total number of sites per formula unit.
"""
total_sites(phase::Phase) = sum(phase.sites)

# =============================================================================
# DSL Functions for defining parameters
# =============================================================================

# Import ForwardDiff for AD compatibility check
using ForwardDiff

"""
    validate_ad_compatibility(func::Function; test_point::Real=1000.0) -> Bool

Check if a function is compatible with ForwardDiff automatic differentiation.

# Arguments

  - `func`: A function `T -> value` to test
  - `test_point`: Temperature value to test at (default: 1000.0 K)

# Returns

  - `true` if the function can be differentiated by ForwardDiff
  - `false` otherwise

# Example

```julia
ghser_ag(T) = -7209.512 + 118.202*T - 23.8463*T*log(T)
validate_ad_compatibility(ghser_ag)  # returns true

# Function with type restriction (may fail)
bad_func(T::Float64) = T^2  # Dual type won't match Float64
validate_ad_compatibility(bad_func)  # returns false
```
"""
function validate_ad_compatibility(func::Function; test_point::Real = 1000.0)
    try
        ForwardDiff.derivative(func, Float64(test_point))
        return true
    catch e
        return false
    end
end

"""
    set_G!(phase::Phase, constituents::Vector{Symbol}, func::Function)

Set the Gibbs energy function for an endmember.

# Arguments

  - `phase`: The phase to modify
  - `constituents`: Species on each sublattice (one per sublattice)
  - `func`: Julia function `T -> G(T)` returning Gibbs energy in J/mol

# Example

```julia
set_G!(fcc, [:AG], T -> -7209.512 + 118.202*T - 23.8463*T*log(T))
```

# Note

The function is automatically checked for ForwardDiff compatibility.
A warning is issued if the function may not work with automatic differentiation.    # Check AD compatibility
"""
function set_G!(phase::Phase, constituents::Vector{Symbol}, func::Function)
    # Check AD compatibility
    if !validate_ad_compatibility(func)
        cons_str = join(string.(constituents), ",")
        @warn "G function for ($cons_str) may not be compatible with automatic differentiation"
    end

    # Wrap each symbol in a vector (one species per sublattice)
    cons = [[c] for c in constituents]
    param = JuliaParameter('G', cons, 0, func)
    push!(phase.julia_parameters, param)
    return phase
end

"""
    set_L!(phase::Phase, constituents::Vector{Symbol}, order::Int, func::Function)

Set a Redlich-Kister interaction parameter.

# Arguments

  - `phase`: The phase to modify
  - `constituents`: Interacting species (e.g., [:AG, :CU] for binary)
  - `order`: Interaction order (0, 1, 2, ...)
  - `func`: Julia function `T -> L(T)` returning interaction energy in J/mol

# Example

```julia
set_L!(fcc, [:AG, :CU], 0, T -> -36583.0 + 4.46*T)
set_L!(fcc, [:AG, :CU], 1, T -> -5765.0 - 4.03*T)
```

# Multi-component systems

For ternary and higher order systems, use the same syntax:

```julia
set_L!(fcc, [:AL, :CU, :MG], 0, T -> -30000.0)  # Ternary interaction
```

# Note    # Check AD compatibility

The function is automatically checked for ForwardDiff compatibility.
"""
function set_L!(phase::Phase, constituents::Vector{Symbol}, order::Int, func::Function)
    # Check AD compatibility
    if !validate_ad_compatibility(func)
        cons_str = join(string.(constituents), ",")
        @warn "L function for ($cons_str;$order) may not be compatible with automatic differentiation"
    end

    # For single sublattice: interaction is within that sublattice
    cons = [constituents]
    param = JuliaParameter('L', cons, order, func)
    push!(phase.julia_parameters, param)
    return phase
end

"""
    set_L!(phase::Phase, constituents::Vector{Vector{Symbol}}, order::Int, func::Function)

Set a Redlich-Kister interaction parameter for multi-sublattice phases.

# Arguments

  - `phase`: The phase to modify
  - `constituents`: Species per sublattice, with interaction in one sublattice
  - `order`: Interaction order (0, 1, 2, ...)
  - `func`: Julia function `T -> L(T)`

# Example (for compound with two sublattices)

```julia
# L(AL,NI:NI;0) - interaction on sublattice 1, NI fixed on sublattice 2
set_L!(phase, [[:AL, :NI], [:NI]], 0, T -> -10000.0)
```
"""
function set_L!(
    phase::Phase,
    constituents::Vector{Vector{Symbol}},
    order::Int,
    func::Function,
)
    # Check AD compatibility
    if !validate_ad_compatibility(func)
        @warn "L function may not be compatible with automatic differentiation"
    end

    param = JuliaParameter('L', constituents, order, func)
    push!(phase.julia_parameters, param)
    return phase
end

"""
    get_julia_endmembers(phase::Phase) -> Vector{JuliaParameter}

Get all endmember (G) parameters defined as Julia functions.
"""
function get_julia_endmembers(phase::Phase)
    return filter(p -> p.symbol == 'G' && is_endmember(p), phase.julia_parameters)
end

"""
    get_julia_interactions(phase::Phase) -> Vector{JuliaParameter}

Get all interaction (L) parameters defined as Julia functions.
"""
function get_julia_interactions(phase::Phase)
    return filter(p -> p.symbol == 'L' || is_interaction(p), phase.julia_parameters)
end

# =============================================================================
# Multi-component support helpers
# =============================================================================

"""
    find_va_index(constituents::Vector{String}) -> Int

Find the index of VA (vacancy) in a sublattice constituent list.

Returns the 1-based index of "VA" in the constituents vector.
If VA is not found, returns 1 as a fallback (for single-constituent sublattices).

# Arguments

  - `constituents`: Vector of species names in a sublattice (e.g., ["C", "VA"])

# Returns

  - Index of "VA" in the vector, or 1 if not found

# Examples

```julia
find_va_index(["VA"])         # => 1
find_va_index(["C", "VA"])    # => 2  (steel1.TDB BCC_A2 sublattice 2)
find_va_index(["CR", "FE"])   # => 1  (fallback, no VA)
```
"""
function find_va_index(constituents::Vector{String})
    idx = findfirst(==("VA"), constituents)
    return isnothing(idx) ? 1 : idx
end

"""
    make_y_matrix(phase::Phase, x::Real, component_indices::Tuple{Int,Int}=(1,2)) -> Matrix

Create a site fraction matrix for a binary calculation within a multi-component phase.

This function properly handles:

 1. Setting the specified two components in the first sublattice
 2. Setting VA (vacancy) to 1.0 in other sublattices, regardless of VA's position

# Arguments

  - `phase`: Phase object with constituent information
  - `x`: Mole fraction of second component (component_indices[2])
  - `component_indices`: Tuple of (first_component_idx, second_component_idx) in sublattice 1
    Default is (1, 2) for standard binary TDB files.

# Returns

  - `y`: Site fraction matrix (n_sublattices × n_max_constituents)

# Examples

```julia
# AgCu FCC_A1 (standard binary)
y = make_y_matrix(fcc, 0.3)  # y[1,1]=0.7 (Ag), y[1,2]=0.3 (Cu), y[2,1]=1.0 (VA)

# Cr-Mo BCC_A2 from steel1.TDB (CR=1, MO=3)
y = make_y_matrix(bcc, 0.5, (1, 3))  # y[1,1]=0.5 (Cr), y[1,3]=0.5 (Mo), y[2,2]=1.0 (VA)
```
"""
function make_y_matrix(
    phase::Phase,
    x::T,
    component_indices::Tuple{Int, Int} = (1, 2),
) where {T <: Real}
    n_sub = length(phase.sites)
    n_const = maximum(length.(phase.constituents))
    y = zeros(T, n_sub, n_const)

    # First sublattice: set specified components
    i1, i2 = component_indices
    y[1, i1] = one(T) - x
    y[1, i2] = x

    # Other sublattices: set VA to 1.0
    for s in 2:n_sub
        va_idx = find_va_index(phase.constituents[s])
        y[s, va_idx] = one(T)
    end

    return y
end
