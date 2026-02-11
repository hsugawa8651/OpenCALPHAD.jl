# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Grid minimizer for composition space scanning

"""
    scan_composition(phase::Phase, T::Real, db::Database, solver::GridSearchSolver;
                     P::Real=1e5, component_indices::Tuple{Int,Int}=(1,2)) -> GridScanResult

Scan Gibbs energy over composition space for a binary phase.

# Arguments

  - `phase`: Phase to scan
  - `T`: Temperature [K]
  - `db`: Database
  - `solver`: GridSearchSolver configuration
  - `P`: Pressure [Pa]
  - `component_indices`: Tuple of (first, second) component indices in sublattice 1.
    Default is (1, 2) for standard binary TDB files.
    For multi-component TDB files (e.g., steel1.TDB), specify the desired components:
    e.g., (1, 3) for Cr-Mo where Cr=1, Mo=3 in BCC_A2.

# Example

```julia
# Standard binary (AgCu)
scan = scan_composition(fcc, 1000.0, db, solver)

# Multi-component (Cr-Mo from steel1.TDB)
scan = scan_composition(bcc, 1800.0, db, solver; component_indices = (1, 3))
```
"""
function scan_composition(
    phase::Phase,
    T::Real,
    db::Database,
    solver::GridSearchSolver;
    P::Real = 1e5,
    component_indices::Tuple{Int, Int} = (1, 2),
)
    n_points = solver.n_points
    min_frac = solver.min_fraction

    # For binary on first sublattice
    x_grid = range(min_frac, 1.0 - min_frac, length = n_points)
    G_values = Float64[]
    y_values = Matrix{Float64}[]

    for x in x_grid
        # Construct site fraction matrix using make_y_matrix
        # Note: make_y_matrix sets y[1, i1] = 1-x and y[1, i2] = x
        # So we pass x directly (fraction of second component)
        y = make_y_matrix(phase, x, component_indices)

        # Clamp and normalize
        clamp_site_fractions!(y; min_val = min_frac)

        # Calculate G
        G = calculate_gibbs_energy(phase, T, y, db; P = P)
        push!(G_values, G)
        push!(y_values, copy(y))
    end

    # Find minimum
    min_idx = argmin(G_values)

    return GridScanResult(phase, collect(x_grid), G_values, y_values, min_idx)
end

"""
    scan_composition_multi(phases::Vector{Phase}, T::Real, db::Database,
                           solver::GridSearchSolver; P::Real=1e5)
                           -> Dict{String, GridScanResult}

Scan composition for multiple phases.
"""
function scan_composition_multi(
    phases::Vector{Phase},
    T::Real,
    db::Database,
    solver::GridSearchSolver;
    P::Real = 1e5,
)
    results = Dict{String, GridScanResult}()

    for phase in phases
        # Skip phases that can't be scanned (e.g., no mixing)
        if isempty(phase.constituents) || length(phase.constituents[1]) < 2
            continue
        end

        result = scan_composition(phase, T, db, solver; P = P)
        results[phase.name] = result
    end

    return results
end

"""
    find_minimum_envelope(scan_results::Dict{String, GridScanResult})
                          -> Tuple{Vector{Float64}, Vector{String}, Vector{Float64}}

Find the minimum G envelope across all phases.

Returns:

  - x_grid: Composition grid
  - stable_phases: Phase name at each composition
  - G_min: Minimum G at each composition
"""
function find_minimum_envelope(scan_results::Dict{String, GridScanResult})
    if isempty(scan_results)
        return Float64[], String[], Float64[]
    end

    # Get common x_grid (assume all have same grid)
    first_result = first(values(scan_results))
    x_grid = first_result.x_grid
    n_points = length(x_grid)

    stable_phases = fill("", n_points)
    G_min = fill(Inf, n_points)

    for (name, result) in scan_results
        for (i, G) in enumerate(result.G_values)
            if G < G_min[i]
                G_min[i] = G
                stable_phases[i] = name
            end
        end
    end

    return x_grid, stable_phases, G_min
end

"""
    find_two_phase_region(scan_results::Dict{String, GridScanResult};
                          tol::Float64=100.0) -> Vector{Tuple{Float64, String, String}}

Find regions where two phases have similar G (potential two-phase regions).

Returns list of (composition, phase1, phase2) tuples where |G1 - G2| < tol.
"""
function find_two_phase_region(
    scan_results::Dict{String, GridScanResult};
    tol::Float64 = 100.0,
)
    if length(scan_results) < 2
        return Tuple{Float64, String, String}[]
    end

    regions = Tuple{Float64, String, String}[]

    # Get phase names and results
    phase_names = collect(keys(scan_results))
    results = collect(values(scan_results))

    # Assume same x_grid
    x_grid = results[1].x_grid
    n_points = length(x_grid)

    for i in 1:n_points
        G_values = [(name, scan_results[name].G_values[i]) for name in phase_names]
        sort!(G_values, by = x -> x[2])

        # Check if top two are within tolerance
        if length(G_values) >= 2
            name1, G1 = G_values[1]
            name2, G2 = G_values[2]
            if abs(G2 - G1) < tol
                push!(regions, (x_grid[i], name1, name2))
            end
        end
    end

    return regions
end

"""
    initial_guess_from_grid(phase::Phase, T::Real, db::Database,
                            solver::GridSearchSolver; P::Real=1e5) -> Matrix{Float64}

Get initial site fractions from grid minimum.
"""
function initial_guess_from_grid(
    phase::Phase,
    T::Real,
    db::Database,
    solver::GridSearchSolver;
    P::Real = 1e5,
)
    result = scan_composition(phase, T, db, solver; P = P)
    return result.y_values[result.min_idx]
end

# Legacy types for compatibility
struct GridPoint
    phase_idx::Int
    site_fractions::Vector{Float64}
    gibbs_energy::Float64
end
