# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Database type definition

"""
    Database

Container for thermodynamic database information.

# Fields

  - `elements::Vector{Element}`: Chemical elements
  - `functions::Vector{GFunction}`: Gibbs energy functions
  - `phases::Vector{Phase}`: Phase definitions with parameters
  - `source::String`: Source file path or identifier
"""
struct Database
    elements::Vector{Element}
    functions::Vector{GFunction}
    phases::Vector{Phase}
    source::String
end

function Database()
    return Database(Element[], GFunction[], Phase[], "")
end

function Database(elements::Vector{Element},
    functions::Vector{GFunction},
    phases::Vector{Phase})
    return Database(elements, functions, phases, "")
end

function Base.show(io::IO, db::Database)
    println(io, "Database: $(db.source)")
    println(io, "  Elements: $(length(db.elements))")
    println(io, "  Functions: $(length(db.functions))")
    println(io, "  Phases: $(length(db.phases))")
end

"""
    get_element(db::Database, name::AbstractString) -> Union{Element, Nothing}

Find an element by name in the database.
"""
function get_element(db::Database, name::AbstractString)
    name_upper = uppercase(name)
    idx = findfirst(e -> uppercase(e.name) == name_upper, db.elements)
    return isnothing(idx) ? nothing : db.elements[idx]
end

"""
    get_phase(db::Database, name::AbstractString) -> Union{Phase, Nothing}

Find a phase by name in the database.
"""
function get_phase(db::Database, name::AbstractString)
    name_upper = uppercase(name)
    idx = findfirst(p -> uppercase(p.name) == name_upper, db.phases)
    return isnothing(idx) ? nothing : db.phases[idx]
end

"""
    get_function(db::Database, name::AbstractString) -> Union{GFunction, Nothing}

Find a function by name in the database.
"""
function get_function(db::Database, name::AbstractString)
    name_upper = uppercase(name)
    idx = findfirst(f -> uppercase(f.name) == name_upper, db.functions)
    return isnothing(idx) ? nothing : db.functions[idx]
end

# =============================================================================
# DSL Helper Functions
# =============================================================================

"""
    add_element!(db::Database, name::Symbol; mass::Float64=0.0, ref_state::String="FCC_A1")

Add an element to the database.

# Example

```julia
add_element!(db, :AG, mass = 107.868)
add_element!(db, :CU, mass = 63.546)
```
"""
function add_element!(
    db::Database,
    name::Symbol;
    mass::Float64 = 0.0,
    ref_state::String = "FCC_A1",
)
    elem = Element(String(name), ref_state, mass)
    push!(db.elements, elem)
    return db
end

"""
    add_elements!(db::Database, names::Vector{Symbol})

Add multiple elements to the database.

# Example

```julia
add_elements!(db, [:AG, :CU])
add_elements!(db, [:AL, :CU, :MG])  # Ternary system
```
"""
function add_elements!(db::Database, names::Vector{Symbol})
    for name in names
        add_element!(db, name)
    end
    return db
end

"""
    add_phase!(db::Database, phase::Phase)

Add a phase to the database.

# Example

```julia
fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])
set_G!(fcc, [:AG], ghser_ag)
set_G!(fcc, [:CU], ghser_cu)
set_L!(fcc, [:AG, :CU], 0, T -> -36583.0 + 4.46*T)
add_phase!(db, fcc)
```
"""
function add_phase!(db::Database, phase::Phase)
    push!(db.phases, phase)
    return db
end
