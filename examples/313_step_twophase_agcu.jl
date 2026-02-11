# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: step_temperature with two-phase detection
#
# Demonstrates step_temperature detecting miscibility gap in Ag-Cu FCC.
# At T=1000K, x(Ag)=0.5 is inside the two-phase region.

using OpenCALPHAD
using Plots

println("=" ^ 70)
println("step_temperature with Two-Phase Detection: Ag-Cu FCC")
println("=" ^ 70)

# Load Ag-Cu thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

# Composition in two-phase region
x_Ag = 0.5

# Temperature range (crossing critical temperature ~1350K for this TDB)
T_start = 1000.0  # K (well below Tc)
T_stop = 1500.0   # K (above Tc)
T_step = 25.0     # K

println("\nPhase: FCC_A1")
println("Composition: x(Ag) = $(x_Ag)")
println("Temperature range: $(T_start)K to $(T_stop)K")
println()

# Run step_temperature
println("Running step_temperature...")
result = step_temperature(fcc, db, x_Ag, T_start, T_stop, T_step)

println("  Number of points: $(length(result.points))")
println()

# Print results
println("Results:")
println("  T [K]  |  G [kJ/mol]  |  n_phases  |  Phase region")
println("  -------|--------------|------------|---------------")

for p in result.points
    region = p.n_phases == 1 ? "Single phase" : "Two-phase (miscibility gap)"
    println(
        "  $(Int(p.axis_value))    |  $(round(p.gibbs_energy/1000, digits=2))       |     $(p.n_phases)      |  $(region)",
    )
end

# Find critical temperature (transition from 2-phase to 1-phase)
idx_critical = findfirst(p -> p.n_phases == 1, result.points)
if !isnothing(idx_critical) && idx_critical > 1
    T_critical = result.points[idx_critical].axis_value
    println("\nCritical temperature (approximate): $(T_critical)K")
    println("(For this TDB, Tc is approximately 1350K)")
end

# Create plots
T_values = [p.axis_value for p in result.points]
G_values = [p.gibbs_energy for p in result.points]
n_phases = [p.n_phases for p in result.points]

# Plot 1: G vs T with phase regions
p1 = plot(T_values, G_values ./ 1000,
    xlabel = "Temperature [K]",
    ylabel = "Gibbs Energy [kJ/mol]",
    title = "Ag-Cu FCC at x(Ag)=$(x_Ag)",
    linewidth = 2,
    color = :blue,
    marker = :circle,
    markersize = 4,
    legend = :topright,
    label = "G(T)",
    size = (700, 500),
)

# Highlight two-phase region
two_phase_T = T_values[n_phases .== 2]
two_phase_G = G_values[n_phases .== 2]
if !isempty(two_phase_T)
    scatter!(p1, two_phase_T, two_phase_G ./ 1000,
        color = :red, markersize = 6, label = "Two-phase")
end

# Plot 2: n_phases vs T
p2 = plot(T_values, n_phases,
    xlabel = "Temperature [K]",
    ylabel = "Number of phases",
    title = "Phase count vs Temperature",
    linewidth = 2,
    color = :green,
    marker = :square,
    markersize = 4,
    ylims = (0.5, 2.5),
    yticks = [1, 2],
    legend = false,
    size = (700, 300),
)

hline!(p2, [1.5], linestyle = :dash, color = :gray)

# Combine
p_combined = plot(p1, p2,
    layout = @layout([a; b{0.3h}]),
    size = (700, 700),
)

# Save
output_path = joinpath(@__DIR__, "313_step_twophase_agcu.png")
savefig(p_combined, output_path)
println("\nSaved to: $output_path")

display(p_combined)

println("\n" * "=" ^ 70)
println("Summary")
println("=" ^ 70)
println("""

This example demonstrates that step_temperature automatically detects
two-phase regions (miscibility gaps) in the calculation.

Key observations:
- At low T (< Tc), n_phases = 2 (inside miscibility gap)
- At high T (> Tc), n_phases = 1 (single phase, solid solution)
- The critical temperature for this TDB is approximately 1350K

The Gibbs energy reported for two-phase regions is the equilibrium
value calculated by find_miscibility_gap (common tangent construction).
""")
