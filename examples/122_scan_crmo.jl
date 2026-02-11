# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Cr-Mo scan_composition with component_indices
#
# Demonstrates the high-level API for multi-component databases.
# Uses component_indices to select binary components from steel1.TDB (6-element).

using OpenCALPHAD
using Plots

# Load steel1 thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "steel1.TDB")
db = read_tdb(tdb_path)

# Get BCC_A2 phase
bcc = get_phase(db, "BCC_A2")

println("=== scan_composition with component_indices ===")
println()

# Show phase structure to determine component indices
println("BCC_A2 sublattice 1: $(bcc.constituents[1])")
println("  CR=1, FE=2, MO=3, SI=4, V=5")
println()

# Temperature
T = 2000.0  # K

# Create solver
solver = GridSearchSolver(n_points = 51)

# Scan composition for Cr-Mo (component_indices=(1,3))
# x represents the mole fraction of the SECOND component (MO)
println("Scanning Cr-Mo composition at T=$(T)K...")
scan = scan_composition(bcc, T, db, solver; component_indices = (1, 3))

println("  Grid points: $(length(scan.x_grid))")
println("  x(Mo) range: $(scan.x_grid[1]) to $(scan.x_grid[end])")
println(
    "  G range: $(round(minimum(scan.G_values)/1000, digits=2)) to $(round(maximum(scan.G_values)/1000, digits=2)) kJ/mol",
)
println()

# Print sample values
println("Sample Gibbs energies:")
println("  x(Mo)  |  G [kJ/mol]")
println("  -------|-------------")
for i in [1, 26, 51]
    println(
        "  $(round(scan.x_grid[i], digits=2))    |  $(round(scan.G_values[i]/1000, digits=2))",
    )
end
println()

# Create plot
p = plot(
    scan.x_grid,
    scan.G_values ./ 1000,
    xlabel = "x(Mo)",
    ylabel = "Gibbs Energy [kJ/mol]",
    title = "Cr-Mo BCC_A2 at T=$(Int(T))K",
    label = "BCC_A2",
    linewidth = 2,
    color = :blue,
    legend = :topright,
    size = (700, 500),
)

# Add minimum point marker
G_min, idx_min = findmin(scan.G_values)
x_min = scan.x_grid[idx_min]
scatter!(p, [x_min], [G_min / 1000], color = :red, markersize = 8, label = "Minimum")
annotate!(p, x_min + 0.1, G_min / 1000 - 0.5,
    text("x=$(round(x_min, digits=2))", 9, :left))

# Save plot
output_path = joinpath(@__DIR__, "122_scan_crmo.png")
savefig(p, output_path)
println("Saved to $output_path")

display(p)

# Comparison: scan all 3 relevant phases
println()
println("=== Comparison of BCC, FCC, HCP for Cr-Mo ===")
println()

fcc = get_phase(db, "FCC_A1")
hcp = get_phase(db, "HCP_A3")

scan_fcc = scan_composition(fcc, T, db, solver; component_indices = (1, 3))
scan_hcp = scan_composition(hcp, T, db, solver; component_indices = (1, 3))

# Plot all phases
p2 = plot(
    xlabel = "x(Mo)",
    ylabel = "Gibbs Energy [kJ/mol]",
    title = "Cr-Mo Phase Comparison at T=$(Int(T))K",
    legend = :topright,
    size = (700, 500),
)

plot!(p2, scan.x_grid, scan.G_values ./ 1000,
    linewidth = 2, color = :blue, label = "BCC_A2")
plot!(p2, scan_fcc.x_grid, scan_fcc.G_values ./ 1000,
    linewidth = 2, color = :green, label = "FCC_A1")
plot!(p2, scan_hcp.x_grid, scan_hcp.G_values ./ 1000,
    linewidth = 2, color = :orange, label = "HCP_A3")

# Save comparison plot
output_path2 = joinpath(@__DIR__, "122_scan_crmo_comparison.png")
savefig(p2, output_path2)
println("Saved to $output_path2")

display(p2)
