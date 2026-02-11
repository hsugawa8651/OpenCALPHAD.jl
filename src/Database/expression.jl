# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - TDB expression parsing and evaluation

"""
    parse_tdb_expression(raw_expr::AbstractString, temp_low::Real)
        -> (temp_breaks::Vector{Float64}, expressions::Vector{String})

Parse a TDB expression string that may contain multiple temperature segments.

Input format: "expr1; T2 Y expr2; T3 N" or "expr1; T2 Y"

  - Segments separated by "; <temperature> [YN]"
  - Y = more segments follow, N = last segment

Returns:

  - temp_breaks: Temperature breakpoints [T_low, T2, T3, ...]
  - expressions: Expression strings for each segment
"""
function parse_tdb_expression(raw_expr::AbstractString, temp_low::Real)
    temp_breaks = Float64[temp_low]
    expressions = String[]

    # Pattern to match: expression followed by ; temperature Y/N
    # Example: "-8856+157*T; 2180 Y -34869+344*T; 6000 N"
    segment_pattern = r";\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)\s*([YN])"i

    remaining = strip(raw_expr)
    last_pos = 1

    for m in eachmatch(segment_pattern, remaining)
        # Extract expression before this marker
        expr = strip(remaining[last_pos:(m.offset - 1)])
        push!(expressions, expr)

        # Extract temperature
        temp = parse(Float64, m.captures[1])
        push!(temp_breaks, temp)

        # Update position
        last_pos = m.offset + length(m.match)
    end

    # Handle remaining expression (if no markers found, use entire string)
    if isempty(expressions)
        push!(expressions, strip(remaining))
        # Default high temperature if none specified
        push!(temp_breaks, 6000.0)
    end

    return (temp_breaks, expressions)
end

"""
    convert_tdb_to_julia(expr::AbstractString) -> String

Convert TDB expression syntax to Julia syntax.

Conversions:

  - T**2 -> T^2
  - LN(x) -> log(x)
  - EXP(x) -> exp(x)
  - LOG10(x) -> log10(x)
"""
function convert_tdb_to_julia(expr::AbstractString)
    result = expr

    # Power operator: ** -> ^
    result = replace(result, r"\*\*" => "^")

    # Natural log: LN(...) -> log(...)
    result = replace(result, r"\bLN\s*\("i => "log(")

    # Exponential: EXP(...) -> exp(...)
    result = replace(result, r"\bEXP\s*\("i => "exp(")

    # Log base 10: LOG10(...) -> log10(...)
    result = replace(result, r"\bLOG10\s*\("i => "log10(")

    # Handle function references: +GHSERCR# -> GHSERCR
    # (These will be resolved later during evaluation)
    result = replace(result, r"#" => "")

    return result
end

"""
    compile_expression(expr::AbstractString) -> Function

Compile a Julia expression string into a function (T, P) -> value.
"""
function compile_expression(expr::AbstractString)
    julia_expr = convert_tdb_to_julia(expr)

    # Create function that takes T and P and returns the evaluated expression
    # Use invokelatest to handle world age issues with eval
    try
        parsed = Meta.parse("(T, P) -> " * julia_expr)
        f = eval(parsed)
        return (T, P) -> Base.invokelatest(f, T, P)
    catch e
        @warn "Failed to compile expression: $expr" exception=e
        return (T, P) -> NaN
    end
end

"""
    compile_expression_with_functions(expr::AbstractString, db::Database) -> Function

Compile a TDB expression that may contain function references (e.g., GHSERCR).
Returns a function (T, P) -> value that evaluates function references at runtime.
"""
function compile_expression_with_functions(expr::AbstractString, db::Database)
    julia_expr = convert_tdb_to_julia(expr)

    # Find all function references in the expression
    # Pattern: word characters that look like function names (all caps, often starting with G)
    func_refs = Set{String}()
    for m in eachmatch(r"\b([A-Z][A-Z0-9_]+)\b", julia_expr)
        name = m.captures[1]
        # Check if this is a known function in the database
        if !isnothing(get_function(db, name))
            push!(func_refs, name)
        end
    end

    if isempty(func_refs)
        # No function references, use simple compilation
        return compile_expression(expr)
    end

    # Compile all referenced functions (with recursive resolution)
    compiled_funcs = Dict{String, Function}()
    for fname in func_refs
        gfunc = get_function(db, fname)
        if !isnothing(gfunc)
            compile_gfunction_with_db!(gfunc, db)
            compiled_funcs[fname] = gfunc.compiled
        end
    end

    # Build expression with function calls
    modified_expr = julia_expr
    for fname in func_refs
        # Replace function name with a call to the compiled function
        modified_expr =
            replace(modified_expr, Regex("\\b$fname\\b") => "funcs[\"$fname\"](T, P)")
    end

    # Compile the modified expression
    try
        parsed = Meta.parse("(T, P, funcs) -> " * modified_expr)
        f = eval(parsed)
        # Return a closure that includes the function dictionary
        # Use invokelatest to handle world age issues
        return (T, P) -> Base.invokelatest(f, T, P, compiled_funcs)
    catch e
        @warn "Failed to compile expression with functions: $expr" exception=e
        return (T, P) -> NaN
    end
end

"""
    compile_gfunction_with_db!(gfunc::GFunction, db::Database)

Compile the GFunction expressions with function reference resolution.
"""
function compile_gfunction_with_db!(gfunc::GFunction, db::Database)
    if !isnothing(gfunc.compiled)
        return  # Already compiled
    end

    n_segments = length(gfunc.expressions)

    if n_segments == 0
        gfunc.compiled = (T, P) -> 0.0
        return
    end

    if n_segments == 1
        gfunc.compiled = compile_expression_with_functions(gfunc.expressions[1], db)
        return
    end

    # Multiple segments - need to select based on temperature
    compiled_segments =
        [compile_expression_with_functions(expr, db) for expr in gfunc.expressions]
    breaks = gfunc.temp_breaks

    gfunc.compiled = function (T, P)
        for i in 1:length(compiled_segments)
            if i == length(compiled_segments) || T < breaks[i + 1]
                return compiled_segments[i](T, P)
            end
        end
        return compiled_segments[end](T, P)
    end
end

"""
    evaluate_gfunction(gfunc::GFunction, T::Real, P::Real=1e5) -> Float64

Evaluate a GFunction at temperature T and pressure P, selecting the appropriate
temperature segment.
"""
function evaluate_gfunction(gfunc::GFunction, T::Real, P::Real = 1e5)
    # Compile if not already done
    if isnothing(gfunc.compiled)
        compile_gfunction!(gfunc)
    end

    # Call the compiled function
    return gfunc.compiled(T, P)
end

"""
    compile_gfunction!(gfunc::GFunction)

Compile the GFunction expressions into an evaluatable function (T, P) -> value.
"""
function compile_gfunction!(gfunc::GFunction)
    n_segments = length(gfunc.expressions)

    if n_segments == 0
        gfunc.compiled = (T, P) -> 0.0
        return
    end

    if n_segments == 1
        # Single segment - simple case
        gfunc.compiled = compile_expression(gfunc.expressions[1])
        return
    end

    # Multiple segments - need to select based on temperature
    compiled_segments = [compile_expression(expr) for expr in gfunc.expressions]
    breaks = gfunc.temp_breaks

    gfunc.compiled = function (T, P)
        # Find the right segment
        for i in 1:length(compiled_segments)
            if i == length(compiled_segments) || T < breaks[i + 1]
                return compiled_segments[i](T, P)
            end
        end
        return compiled_segments[end](T, P)
    end
end

"""
    evaluate_parameter(param::Parameter, T::Real, db::Database, P::Real=1e5) -> Float64

Evaluate a Parameter at temperature T and pressure P, resolving function references.
"""
function evaluate_parameter(param::Parameter, T::Real, db::Database, P::Real = 1e5)
    n_segments = length(param.expressions)

    if n_segments == 0
        return 0.0
    end

    # Find the right temperature segment
    segment_idx = 1
    for i in 1:(length(param.temp_breaks) - 1)
        if T >= param.temp_breaks[i] && T < param.temp_breaks[i + 1]
            segment_idx = i
            break
        end
        if i == length(param.temp_breaks) - 1
            segment_idx = min(i, n_segments)
        end
    end

    expr = param.expressions[min(segment_idx, n_segments)]
    compiled = compile_expression_with_functions(expr, db)
    return compiled(T, P)
end
