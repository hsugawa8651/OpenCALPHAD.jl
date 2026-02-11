# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Gibbs energy function definitions

"""
    AbstractGFunction

Abstract type for Gibbs energy functions.
"""
abstract type AbstractGFunction end

"""
    GFunction

Represents a FUNCTION definition in TDB format.
Temperature-dependent Gibbs energy expression.

# Fields

  - `name::String`: Function identifier (e.g., "GHSERAG")
  - `temp_range::Tuple{Float64,Float64}`: Valid temperature range [K]
  - `temp_breaks::Vector{Float64}`: Temperature breakpoints for piecewise functions
  - `expressions::Vector{String}`: Expression strings for each temperature segment
  - `compiled::Union{Nothing,Function}`: Compiled Julia function (lazy evaluation)
"""
mutable struct GFunction <: AbstractGFunction
    name::String
    temp_range::Tuple{Float64, Float64}
    temp_breaks::Vector{Float64}
    expressions::Vector{String}
    compiled::Union{Nothing, Function}
end

function GFunction(name::AbstractString,
    temp_range::Tuple{<:Real, <:Real},
    expressions::Vector{<:AbstractString})
    return GFunction(
        String(name),
        (Float64(temp_range[1]), Float64(temp_range[2])),
        Float64[],
        String.(expressions),
        nothing,
    )
end

function GFunction(name::AbstractString,
    temp_breaks::Vector{<:Real},
    expressions::Vector{<:AbstractString})
    temp_range = (minimum(temp_breaks), maximum(temp_breaks))
    return GFunction(
        String(name),
        temp_range,
        Float64.(temp_breaks),
        String.(expressions),
        nothing,
    )
end

function Base.show(io::IO, f::GFunction)
    print(io, "GFunction($(f.name), T=$(f.temp_range[1])-$(f.temp_range[2]))")
end

"""
    evaluate(f::GFunction, T::Real) -> Float64

Evaluate the Gibbs energy function at temperature T.
"""
function evaluate(f::GFunction, T::Real)
    # TODO: Implement proper evaluation with compiled function
    error("GFunction evaluation not yet implemented")
end
