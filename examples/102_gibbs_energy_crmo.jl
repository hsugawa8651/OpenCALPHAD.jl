# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Cr-Mo Gibbs energy curves
#
# Plot Gibbs energy curves for Cr-Mo binary system using steel1.TDB
#
# This example demonstrates two approaches to handle multi-component databases:
# 1. Manual y-matrix construction (shown below) for understanding the internals
# 2. High-level API with component_indices (see examples/122_scan_crmo.jl)
#
# For production use, the high-level API is recommended:
#   scan_composition(bcc, T, db, solver; component_indices=(1, 3))  # CR=1, MO=3

using OpenCALPHAD
using Plots

# Load steel1 thermodynamic database (6-component: C, Cr, Fe, Mo, Si, V)
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "steel1.TDB")
db = read_tdb(tdb_path)

# Get phases
bcc = get_phase(db, "BCC_A2")
liquid = get_phase(db, "LIQUID")

println("=== Cr-Mo System Analysis ===")
println()

# Show phase structures
println("BCC_A2 structure:")
println("  Sublattice 1: $(bcc.constituents[1])")
println("  Sublattice 2: $(bcc.constituents[2])")
println("  Sites: $(bcc.sites)")
println()

println("LIQUID structure:")
println("  Sublattice 1: $(liquid.constituents[1])")
println()

# Component indices (determined from constituent arrays)
# BCC_A2 sublattice 1: CR(1), FE(2), MO(3), SI(4), V(5)
# BCC_A2 sublattice 2: C(1), VA(2)
# LIQUID sublattice 1: C(1), CR(2), FE(3), MO(4), SI(5), V(6)

const CR_IDX_BCC = 1
const MO_IDX_BCC = 3
const VA_IDX_BCC = 2  # in sublattice 2

const CR_IDX_LIQ = 2
const MO_IDX_LIQ = 4

# Temperature for calculation
T = 2000.0  # K (both Cr and Mo are solid at this T, but close to Cr melting point 2180K)

println("Calculating Gibbs energy curves at T=$(T)K...")
println()

# Composition range
x_Mo_values = 0.0:0.02:1.0
n_points = length(x_Mo_values)

# Calculate G for BCC_A2
G_bcc = zeros(n_points)
n_sub_bcc = length(bcc.sites)
n_const_bcc = length(bcc.constituents[1])

for (i, x_Mo) in enumerate(x_Mo_values)
    y = zeros(n_sub_bcc, max(n_const_bcc, length(bcc.constituents[2])))
    y[1, CR_IDX_BCC] = 1.0 - x_Mo  # Cr
    y[1, MO_IDX_BCC] = x_Mo        # Mo
    y[2, VA_IDX_BCC] = 1.0         # VA (fixed)

    G_bcc[i] = calculate_gibbs_energy(bcc, T, y, db)
end

# Calculate G for LIQUID
G_liquid = zeros(n_points)
n_const_liq = length(liquid.constituents[1])

for (i, x_Mo) in enumerate(x_Mo_values)
    y = zeros(1, n_const_liq)
    y[1, CR_IDX_LIQ] = 1.0 - x_Mo  # Cr
    y[1, MO_IDX_LIQ] = x_Mo        # Mo

    G_liquid[i] = calculate_gibbs_energy(liquid, T, y, db)
end

# Print some values
println("Sample Gibbs energies at T=$(T)K:")
println("  x(Mo)  |  G_BCC [kJ/mol]  |  G_LIQUID [kJ/mol]  |  G_LIQ - G_BCC")
println("  -------|------------------|---------------------|----------------")
for i in [1, div(n_points, 2), n_points]
    x = x_Mo_values[i]
    dG = G_liquid[i] - G_bcc[i]
    println(
        "  $(round(x, digits=2))    |  $(round(G_bcc[i]/1000, digits=2))         |  $(round(G_liquid[i]/1000, digits=2))            |  $(round(dG/1000, digits=2))",
    )
end
println()

# Create plot
p = plot(
    xlabel = "x(Mo)",
    ylabel = "Gibbs Energy [kJ/mol]",
    title = "Cr-Mo Gibbs Energy at T=$(Int(T))K",
    legend = :topright,
    size = (700, 500),
)

# Plot BCC_A2
plot!(p, collect(x_Mo_values), G_bcc ./ 1000,
    linewidth = 2,
    color = :blue,
    label = "BCC_A2",
)

# Plot LIQUID
plot!(p, collect(x_Mo_values), G_liquid ./ 1000,
    linewidth = 2,
    color = :red,
    label = "LIQUID",
)

# Add annotations
annotate!(p, 0.5, minimum(G_bcc)/1000 - 2,
    text("BCC is stable\n(G_BCC < G_LIQUID)", 10, :center))

# Save to file
output_path = joinpath(@__DIR__, "102_gibbs_energy_crmo.png")
savefig(p, output_path)
println("Saved to $output_path")

# Display plot
display(p)

# Additional: Temperature dependence
println()
println("=== Temperature Dependence at x(Mo)=0.5 ===")
println()

T_values = 1500.0:100.0:2500.0
G_bcc_T = zeros(length(T_values))
G_liquid_T = zeros(length(T_values))

x_Mo = 0.5

# BCC y-matrix for x_Mo = 0.5
y_bcc = zeros(n_sub_bcc, max(n_const_bcc, length(bcc.constituents[2])))
y_bcc[1, CR_IDX_BCC] = 1.0 - x_Mo
y_bcc[1, MO_IDX_BCC] = x_Mo
y_bcc[2, VA_IDX_BCC] = 1.0

# LIQUID y-matrix for x_Mo = 0.5
y_liq = zeros(1, n_const_liq)
y_liq[1, CR_IDX_LIQ] = 1.0 - x_Mo
y_liq[1, MO_IDX_LIQ] = x_Mo

for (i, T_i) in enumerate(T_values)
    G_bcc_T[i] = calculate_gibbs_energy(bcc, T_i, y_bcc, db)
    G_liquid_T[i] = calculate_gibbs_energy(liquid, T_i, y_liq, db)
end

# Find approximate melting temperature (where G_BCC ≈ G_LIQUID)
dG = G_liquid_T .- G_bcc_T
idx_cross = findfirst(dG .< 0)
if !isnothing(idx_cross) && idx_cross > 1
    # Linear interpolation
    T1, T2 = T_values[idx_cross - 1], T_values[idx_cross]
    dG1, dG2 = dG[idx_cross - 1], dG[idx_cross]
    T_melt = T1 - dG1 * (T2 - T1) / (dG2 - dG1)
    println("Approximate melting temperature at x(Mo)=0.5: $(round(T_melt, digits=0))K")
else
    println("No crossing found in temperature range")
end

# Plot G vs T
p2 = plot(
    xlabel = "Temperature [K]",
    ylabel = "Gibbs Energy [kJ/mol]",
    title = "Cr-Mo (x_Mo=0.5) G vs T",
    legend = :topright,
    size = (700, 500),
)

plot!(p2, collect(T_values), G_bcc_T ./ 1000,
    linewidth = 2,
    color = :blue,
    label = "BCC_A2",
)

plot!(p2, collect(T_values), G_liquid_T ./ 1000,
    linewidth = 2,
    color = :red,
    label = "LIQUID",
)

# Save G-T plot
output_path2 = joinpath(@__DIR__, "102_gibbs_energy_crmo_T.png")
savefig(p2, output_path2)
println("Saved to $output_path2")

display(p2)
