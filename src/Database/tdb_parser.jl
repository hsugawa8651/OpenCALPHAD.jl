# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - TDB file parser

"""
    read_tdb(filepath::AbstractString) -> Database

Parse a TDB (Thermodynamic Database) file and return a Database object.

# Arguments

  - `filepath`: Path to the TDB file

# Returns

  - `Database`: Parsed database containing elements, functions, phases, and parameters
"""
function read_tdb(filepath::AbstractString)
    if !isfile(filepath)
        error("TDB file not found: $filepath")
    end

    content = read(filepath, String)
    return parse_tdb(content, filepath)
end

"""
    parse_tdb(content::AbstractString, source::AbstractString="") -> Database

Parse TDB content string and return a Database object.
"""
function parse_tdb(content::AbstractString, source::AbstractString = "")
    elements = Element[]
    functions = GFunction[]
    phases = Phase[]

    # Normalize line endings and join continuation lines
    lines = preprocess_tdb(content)

    i = 1
    while i <= length(lines)
        line = strip(lines[i])

        if isempty(line) || startswith(line, "\$")
            # Skip empty lines and comments
            i += 1
            continue
        end

        keyword = get_keyword(line)

        if keyword == "ELEMENT"
            elem = parse_element(line)
            if !isnothing(elem)
                push!(elements, elem)
            end
        elseif keyword == "FUNCTION"
            func, i = parse_function(lines, i)
            if !isnothing(func)
                push!(functions, func)
            end
        elseif keyword == "PHASE"
            phase = parse_phase(line)
            if !isnothing(phase)
                push!(phases, phase)
            end
        elseif keyword == "CONSTITUENT"
            parse_constituent!(phases, line)
        elseif keyword == "PARAMETER"
            param, i = parse_parameter(lines, i)
            if !isnothing(param)
                add_parameter!(phases, param)
            end
        end

        i += 1
    end

    return Database(elements, functions, phases, source)
end

"""
    preprocess_tdb(content::AbstractString) -> Vector{String}

Preprocess TDB content: normalize line endings and join continuation lines.
"""
function preprocess_tdb(content::AbstractString)
    # Normalize line endings
    content = replace(content, "\r\n" => "\n")
    content = replace(content, "\r" => "\n")

    # Split into lines
    raw_lines = split(content, "\n")

    # Join continuation lines (lines ending with content before !)
    lines = String[]
    current_line = ""

    for line in raw_lines
        stripped = strip(line)
        if isempty(stripped)
            if !isempty(current_line)
                push!(lines, current_line)
                current_line = ""
            end
            continue
        end

        current_line *= " " * stripped

        # Check if line is complete (ends with !)
        if endswith(stripped, "!")
            push!(lines, strip(current_line))
            current_line = ""
        end
    end

    if !isempty(current_line)
        push!(lines, strip(current_line))
    end

    return lines
end

"""
    get_keyword(line::AbstractString) -> String

Extract the keyword from a TDB line.
"""
function get_keyword(line::AbstractString)
    m = match(r"^\s*(\w+)", line)
    return isnothing(m) ? "" : uppercase(m.captures[1])
end

"""
    parse_element(line::AbstractString) -> Union{Element, Nothing}

Parse an ELEMENT line from TDB.
Format: ELEMENT <name> <reference> <mass> <H298> <S298> !
"""
function parse_element(line::AbstractString)
    # Remove trailing !
    line = replace(line, r"!\s*$" => "")

    # Match ELEMENT line pattern
    pattern = r"ELEMENT\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)"i
    m = match(pattern, line)

    if isnothing(m)
        return nothing
    end

    name = uppercase(m.captures[1])
    reference = uppercase(m.captures[2])
    mass = parse(Float64, m.captures[3])
    H298 = parse(Float64, m.captures[4])
    S298 = parse(Float64, m.captures[5])

    return Element(name, reference, mass, H298, S298)
end

"""
    parse_phase(line::AbstractString) -> Union{Phase, Nothing}

Parse a PHASE line from TDB.
Format: PHASE <name> <type> <n_sublattices> <site1> <site2> ... !
"""
function parse_phase(line::AbstractString)
    line = replace(line, r"!\s*$" => "")

    # Extract phase name, type, and sites
    pattern = r"PHASE\s+(\S+)\s+(\S+)\s+(\d+)\s+(.*)"i
    m = match(pattern, line)

    if isnothing(m)
        return nothing
    end

    name = uppercase(m.captures[1])
    # Extract state from name (e.g., LIQUID:L -> 'L')
    state = 'S'
    if occursin(":", name)
        parts = split(name, ":")
        name = parts[1]
        if length(parts) > 1 && length(parts[2]) == 1
            state = uppercase(parts[2])[1]
        end
    end
    if occursin("LIQUID", name)
        state = 'L'
    elseif occursin("GAS", name)
        state = 'G'
    end

    model_type = m.captures[2]
    # n_sublattices = parse(Int, m.captures[3])
    sites_str = m.captures[4]

    # Parse site ratios
    sites = Float64[]
    for s in split(sites_str)
        s = strip(s)
        if !isempty(s)
            try
                push!(sites, parse(Float64, s))
            catch
                # Ignore non-numeric entries
            end
        end
    end

    return Phase(name, state, model_type, sites)
end

"""
    parse_constituent!(phases::Vector{Phase}, line::AbstractString)

Parse a CONSTITUENT line and update the corresponding phase.
Format: CONSTITUENT <phase> :<sub1>,<sub2>,...: <species1>,<species2>:... !
"""
function parse_constituent!(phases::Vector{Phase}, line::AbstractString)
    line = replace(line, r"!\s*$" => "")

    pattern = r"CONSTITUENT\s+(\S+)\s*:\s*(.*)"i
    m = match(pattern, line)

    if isnothing(m)
        return
    end

    phase_name = uppercase(split(m.captures[1], ":")[1])
    constituents_str = m.captures[2]

    # Find the phase
    idx = findfirst(p -> uppercase(p.name) == phase_name, phases)
    if isnothing(idx)
        return
    end

    # Parse constituents for each sublattice
    sublattices = split(constituents_str, ":")
    constituents = Constitution()

    for sub_str in sublattices
        sub_str = strip(sub_str)
        if isempty(sub_str)
            continue
        end
        # Strip trailing % (species type indicator) and normalize
        species = [uppercase(replace(strip(s), r"%$" => "")) for s in split(sub_str, ",")]
        push!(constituents, species)
    end

    phases[idx].constituents = constituents
end

"""
    parse_function(lines::Vector{String}, start_idx::Int) -> Tuple{Union{GFunction, Nothing}, Int}

Parse a FUNCTION definition from TDB.
Handles multiple temperature segments: "expr1; T2 Y expr2; T3 N"
"""
function parse_function(lines::Vector{String}, start_idx::Int)
    line = lines[start_idx]
    line = replace(line, r"!\s*$" => "")

    # Pattern: FUNCTION name temp_low expression...
    pattern = r"FUNCTION\s+(\S+)\s+(\S+)\s+(.*)"i
    m = match(pattern, line)

    if isnothing(m)
        return (nothing, start_idx)
    end

    name = uppercase(m.captures[1])
    temp_low = parse(Float64, m.captures[2])
    raw_expr = strip(m.captures[3])

    # Parse temperature segments using the expression parser
    # Note: parse_tdb_expression is defined in expression.jl (loaded after this file)
    # We'll do inline parsing here for now

    temp_breaks = Float64[temp_low]
    expressions = String[]

    # Pattern: expression followed by ; temperature Y/N
    segment_pattern = r";\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)\s*([YN])"i

    remaining = raw_expr
    last_pos = 1

    for m_seg in eachmatch(segment_pattern, remaining)
        # Extract expression before this marker
        expr = strip(remaining[last_pos:(m_seg.offset - 1)])
        push!(expressions, expr)

        # Extract temperature
        temp = parse(Float64, m_seg.captures[1])
        push!(temp_breaks, temp)

        # Update position
        last_pos = m_seg.offset + length(m_seg.match)
    end

    # Handle remaining expression after last marker
    if last_pos <= length(remaining)
        final_expr = strip(remaining[last_pos:end])
        if !isempty(final_expr)
            push!(expressions, final_expr)
        end
    end

    # If no segments found, use entire string as single expression
    if isempty(expressions)
        push!(expressions, strip(remaining))
        push!(temp_breaks, 6000.0)  # Default high temperature
    end

    temp_range = (temp_breaks[1], temp_breaks[end])
    return (GFunction(name, temp_range, temp_breaks, expressions, nothing), start_idx)
end

"""
    parse_parameter(lines::Vector{String}, start_idx::Int) -> Tuple{Union{Parameter, Nothing}, Int}

Parse a PARAMETER definition from TDB.
Handles multiple temperature segments: "expr1; T2 Y expr2; T3 N"
"""
function parse_parameter(lines::Vector{String}, start_idx::Int)
    line = lines[start_idx]
    line = replace(line, r"!\s*$" => "")

    # Remove trailing reference (e.g., "REF283")
    line = replace(line, r"\s+REF\d+\s*$"i => "")

    # Pattern: PARAMETER G/L/TC/BMAGN(PHASE,A:B:C;order) temp expression...
    pattern =
        r"PARAMETER\s+([GLTB]\w*)\s*\(\s*(\S+?)\s*,\s*([^;)]+)\s*(?:;\s*(\d+))?\s*\)\s+(\S+)\s+(.*)"i
    m = match(pattern, line)

    if isnothing(m)
        return (nothing, start_idx)
    end

    symbol = uppercase(m.captures[1])[1]
    phase_name = uppercase(split(m.captures[2], ":")[1])
    constituents_str = m.captures[3]
    order = isnothing(m.captures[4]) ? 0 : parse(Int, m.captures[4])
    temp_low = parse(Float64, m.captures[5])
    raw_expr = strip(m.captures[6])

    # Parse constituents
    sublattices = split(constituents_str, ":")
    constituents = Constitution()
    for sub_str in sublattices
        sub_str = strip(sub_str)
        species = [uppercase(strip(s)) for s in split(sub_str, ",")]
        push!(constituents, species)
    end

    # Parse temperature segments (same logic as parse_function)
    temp_breaks = Float64[temp_low]
    expressions = String[]

    segment_pattern = r";\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)\s*([YN])"i

    remaining = raw_expr
    last_pos = 1

    for m_seg in eachmatch(segment_pattern, remaining)
        expr = strip(remaining[last_pos:(m_seg.offset - 1)])
        push!(expressions, expr)
        temp = parse(Float64, m_seg.captures[1])
        push!(temp_breaks, temp)
        last_pos = m_seg.offset + length(m_seg.match)
    end

    # Handle remaining expression after last marker
    if last_pos <= length(remaining)
        final_expr = strip(remaining[last_pos:end])
        if !isempty(final_expr)
            push!(expressions, final_expr)
        end
    end

    # If no segments found, use entire string as single expression
    if isempty(expressions)
        push!(expressions, strip(remaining))
        push!(temp_breaks, 6000.0)  # Default high temperature
    end

    return (
        Parameter(phase_name, symbol, constituents, order, temp_breaks, expressions),
        start_idx,
    )
end

"""
    add_parameter!(phases::Vector{Phase}, param::Parameter)

Add a parameter to the corresponding phase.
"""
function add_parameter!(phases::Vector{Phase}, param::Parameter)
    idx = findfirst(p -> uppercase(p.name) == uppercase(param.phase_name), phases)
    if !isnothing(idx)
        push!(phases[idx].parameters, param)
    end
end
