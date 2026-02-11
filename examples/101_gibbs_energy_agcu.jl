# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Gibbs energy curve with common tangent
#
# Plot Gibbs energy curve with common tangent construction
#
# This script demonstrates how to calculate and plot the Gibbs energy curve
# for a binary system, showing the common tangent construction that defines
# the miscibility gap.

using OpenCALPHAD
using Plots

# Load Ag-Cu thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

# Temperature for calculation
T = 1000.0  # K

println("Calculating Gibbs energy curve at T=$(T)K...")

# Scan composition to get G(x) curve
solver = GridSearchSolver(n_points = 101)
scan = scan_composition(fcc, T, db, solver)

# Create plot (convert J/mol to kJ/mol)
p = plot(scan.x_grid, scan.G_values ./ 1000,
    xlabel = "x(Ag)",
    ylabel = "Gibbs Energy [kJ/mol]",
    title = "Ag-Cu FCC Gibbs Energy at T=$(Int(T))K",
    linewidth = 2,
    legend = :topright,
    label = "G(x)",
    size = (600, 400),
)

# Find and plot the common tangent (miscibility gap)
gap = find_miscibility_gap(fcc, T, db)
if !isnothing(gap)
    println("Miscibility gap found:")
    println("  x1 = $(round(gap.x1, digits=4)) (Cu-rich)")
    println("  x2 = $(round(gap.x2, digits=4)) (Ag-rich)")

    # Mark equilibrium compositions
    scatter!(p, [gap.x1, gap.x2], [gap.G1, gap.G2] ./ 1000,
        markersize = 8,
        color = :red,
        label = "Equilibrium compositions",
    )

    # Draw common tangent line
    plot!(p, [gap.x1, gap.x2], [gap.G1, gap.G2] ./ 1000,
        linestyle = :dash,
        color = :red,
        linewidth = 2,
        label = "Common tangent",
    )
else
    println("No miscibility gap found at this temperature.")
end

# Save to file
output_path = joinpath(@__DIR__, "101_gibbs_energy_agcu.png")
savefig(p, output_path)
println("Saved to $output_path")

# Display plot (if running interactively)
display(p)
