# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Cr-Mo step_temperature with component_indices
#
# Demonstrates step_temperature for multi-component databases.
# Uses component_indices to select Cr-Mo binary from steel1.TDB.

using OpenCALPHAD
using Plots

# Load steel1 thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "steel1.TDB")
db = read_tdb(tdb_path)

# Get BCC_A2 phase
bcc = get_phase(db, "BCC_A2")

println("=== step_temperature with component_indices ===")
println()

# Component indices for Cr-Mo
# BCC_A2 sublattice 1: CR(1), FE(2), MO(3), SI(4), V(5)
CR_IDX = 1
MO_IDX = 3

# Temperature range (limited for faster execution)
T_start = 1800.0  # K
T_stop = 2200.0   # K
T_step = 100.0    # K

# Composition (Mo mole fraction)
x_Mo = 0.5

println("Phase: BCC_A2")
println("Composition: x(Mo) = $(x_Mo)")
println("Temperature range: $(T_start)K to $(T_stop)K, step $(T_step)K")
println()

# Run step_temperature
println("Running step_temperature...")
result = step_temperature(bcc, db, x_Mo, T_start, T_stop, T_step;
    component_indices = (CR_IDX, MO_IDX))

println("  Number of points: $(length(result.points))")
println("  All converged: $(all(p -> p.converged, result.points))")
println()

# Extract data (axis_value is temperature for step_temperature)
T_values = [p.axis_value for p in result.points]
G_values = [p.gibbs_energy for p in result.points]

# Print sample values
println("Sample results:")
println("  T [K]  |  G [kJ/mol]  |  n_phases")
println("  -------|--------------|----------")
for p in result.points
    println(
        "  $(Int(p.axis_value))   |  $(round(p.gibbs_energy/1000, digits=2))       |  $(p.n_phases)",
    )
end
println()

# Create G vs T plot
p1 = plot(
    T_values,
    G_values ./ 1000,
    xlabel = "Temperature [K]",
    ylabel = "Gibbs Energy [kJ/mol]",
    title = "Cr-Mo BCC_A2 at x(Mo)=$(x_Mo)",
    label = "BCC_A2",
    linewidth = 2,
    color = :blue,
    marker = :circle,
    markersize = 3,
    legend = :topright,
    size = (700, 500),
)

# Save plot
output_path = joinpath(@__DIR__, "322_step_crmo.png")
savefig(p1, output_path)
println("Saved to $output_path")

display(p1)

# Note: For comparison at different compositions, see the commented code below.
# It takes longer to run due to multiple step_temperature calls.
#
# x_values = [0.0, 0.25, 0.5, 0.75, 1.0]
# for x in x_values
#     result_i = step_temperature(bcc, db, x, T_start, T_stop, T_step;
#         component_indices = (CR_IDX, MO_IDX))
#     # ... plot
# end
