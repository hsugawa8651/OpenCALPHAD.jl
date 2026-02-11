# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Binary phase diagram
#
# Plot Ag-Cu binary phase diagram
#
# This script demonstrates how to calculate and plot a binary phase diagram
# showing the FCC miscibility gap in the Ag-Cu system.

using OpenCALPHAD
using Plots

# Load Ag-Cu thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

println("Calculating Ag-Cu phase diagram...")

# Calculate phase diagram from 800K to 1300K
result = map_phase_diagram(fcc, db, 800.0, 1300.0, 25.0)

# Extract converged points
converged = filter(p -> p.converged, result.points)
println("Converged points: $(length(converged)) / $(length(result.points))")

# Extract boundary data
T_vals = [p.temperature for p in converged]
x1_vals = [p.compositions[1] for p in converged]
x2_vals = [p.compositions[2] for p in converged]

# Create phase diagram plot
p = plot(x1_vals, T_vals,
    label = "Cu-rich boundary",
    xlabel = "x(Ag)",
    ylabel = "Temperature [K]",
    title = "Ag-Cu Phase Diagram (FCC Miscibility Gap)",
    linewidth = 2,
    legend = :top,
    size = (600, 400),
)
plot!(p, x2_vals, T_vals, label = "Ag-rich boundary", linewidth = 2)

# Fill the two-phase region
x_fill = vcat(x1_vals, reverse(x2_vals))
T_fill = vcat(T_vals, reverse(T_vals))
plot!(p, x_fill, T_fill,
    fillrange = minimum(T_vals),
    fillalpha = 0.2,
    label = "Two-phase region",
    linewidth = 0,
)

# Save to file
output_path = joinpath(@__DIR__, "411_phase_diagram_agcu.png")
savefig(p, output_path)
println("Saved to $output_path")

# Display plot (if running interactively)
display(p)
