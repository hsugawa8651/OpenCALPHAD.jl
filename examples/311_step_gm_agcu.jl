# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Thermodynamic property curves (step2.OCM equivalent)
#
# Plot G, H_mix, S_mix, and G_mix curves for Ag-Cu FCC at fixed temperature.
# This corresponds to step2.OCM from openCALPHAD (Fortran).
#
# Uses the new thermodynamic property functions:
#   - calculate_mixing_gibbs, calculate_mixing_enthalpy, calculate_mixing_entropy
#   - thermodynamic_properties (efficient one-call version)
#
# Usage:
#   julia --project examples/311_step_gm_agcu.jl

using OpenCALPHAD
using Plots

println("=" ^ 70)
println("Thermodynamic Property Curves: Ag-Cu FCC (step2.OCM equivalent)")
println("=" ^ 70)

# Load Ag-Cu thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

# Temperature for calculation
T = 1000.0  # K
P = 1e5     # Pa

println("\nCalculating at T = $(T) K...")

# =============================================================================
# Calculate properties over composition using new API
# =============================================================================

solver = GridSearchSolver(n_points = 101)
scan = scan_composition(fcc, T, db, solver)

# Calculate mixing properties for each composition point
G_values = scan.G_values
G_mix = Float64[]
H_mix = Float64[]
S_mix = Float64[]

for y in scan.y_values
    props = thermodynamic_properties(fcc, T, y, db; P = P)
    push!(G_mix, props.G_mix)
    push!(H_mix, props.H_mix)
    push!(S_mix, props.S_mix)
end

# x_Ag = scan.x_grid (fraction of first constituent = AG)
x_Ag = scan.x_grid
x_Cu = 1.0 .- x_Ag

# Reference line for Gibbs energy plot
n_sub = length(fcc.sites)
n_const = length(fcc.constituents[1])

y_Ag = zeros(n_sub, n_const)
y_Ag[1, 1] = 1.0
for s in 2:n_sub
    y_Ag[s, 1] = 1.0
end
G_Ag_pure = calculate_gibbs_energy(fcc, T, y_Ag, db; P = P)

y_Cu = zeros(n_sub, n_const)
y_Cu[1, 2] = 1.0
for s in 2:n_sub
    y_Cu[s, 1] = 1.0
end
G_Cu_pure = calculate_gibbs_energy(fcc, T, y_Cu, db; P = P)

G_ref_line = x_Ag .* G_Ag_pure .+ x_Cu .* G_Cu_pure

println("Pure Ag (FCC) G = $(round(G_Ag_pure/1000, digits=2)) kJ/mol")
println("Pure Cu (FCC) G = $(round(G_Cu_pure/1000, digits=2)) kJ/mol")

# =============================================================================
# Print values at x = 0.5
# =============================================================================

println("\nAt x(Ag) = 0.5:")
idx_mid = findfirst(x -> x >= 0.5, x_Ag)
println("  G_mix = $(round(G_mix[idx_mid], digits=1)) J/mol")
println("  H_mix = $(round(H_mix[idx_mid], digits=1)) J/mol")
println("  S_mix = $(round(S_mix[idx_mid], digits=3)) J/(mol·K)")
println("  -T*S_mix = $(round(-T*S_mix[idx_mid], digits=1)) J/mol")

# =============================================================================
# Find miscibility gap
# =============================================================================

gap = find_miscibility_gap(fcc, T, db)
if !isnothing(gap)
    println("\nMiscibility gap:")
    println("  x1 = $(round(gap.x1, digits=4)) (Cu-rich)")
    println("  x2 = $(round(gap.x2, digits=4)) (Ag-rich)")
end

# =============================================================================
# Create 2x2 plot panel
# =============================================================================

println("\nGenerating plots...")

# Panel (a): Gibbs energy with common tangent
pa = plot(x_Ag, G_values ./ 1000,
    xlabel = "x(Ag)", ylabel = "G [kJ/mol]",
    title = "(a) Gibbs Energy",
    label = "G(x)", linewidth = 2, color = :blue,
    legend = :topright)

plot!(pa, x_Ag, G_ref_line ./ 1000,
    label = "Reference", linewidth = 1, linestyle = :dash, color = :gray)

if !isnothing(gap)
    scatter!(pa, [gap.x1, gap.x2], [gap.G1, gap.G2] ./ 1000,
        markersize = 6, color = :red, label = "Equilibrium")
    plot!(pa, [gap.x1, gap.x2], [gap.G1, gap.G2] ./ 1000,
        linewidth = 2, linestyle = :dash, color = :red, label = "")
end

# Panel (b): Mixing Gibbs energy
pb = plot(x_Ag, G_mix,
    xlabel = "x(Ag)", ylabel = "G_mix [J/mol]",
    title = "(b) Mixing Gibbs Energy",
    label = "G_mix", linewidth = 2, color = :blue,
    legend = :bottomright)

hline!(pb, [0], linestyle = :dash, color = :gray, label = "")

if !isnothing(gap)
    vspan!(pb, [gap.x1, gap.x2], alpha = 0.2, color = :red, label = "Two-phase")
end

# Panel (c): Mixing enthalpy
pc = plot(x_Ag, H_mix,
    xlabel = "x(Ag)", ylabel = "H_mix [J/mol]",
    title = "(c) Mixing Enthalpy",
    label = "H_mix", linewidth = 2, color = :red,
    legend = :bottomright)

hline!(pc, [0], linestyle = :dash, color = :gray, label = "")

# Panel (d): Mixing entropy
pd = plot(x_Ag, S_mix,
    xlabel = "x(Ag)", ylabel = "S_mix [J/(mol·K)]",
    title = "(d) Mixing Entropy",
    label = "S_mix", linewidth = 2, color = :green,
    legend = :topright)

# Ideal mixing entropy for comparison
R = 8.314462618
S_ideal = -R .* (x_Ag .* log.(max.(x_Ag, 1e-10)) .+ x_Cu .* log.(max.(x_Cu, 1e-10)))
plot!(pd, x_Ag, S_ideal,
    label = "S_ideal", linewidth = 1, linestyle = :dash, color = :orange)

# Combine panels
p_combined = plot(pa, pb, pc, pd,
    layout = (2, 2), size = (900, 700),
    plot_title = "Ag-Cu FCC at T=$(Int(T))K",
    margin = 5Plots.mm)

# Save
output_path = joinpath(@__DIR__, "311_step_gm_agcu.png")
savefig(p_combined, output_path)
println("\nSaved to: $output_path")

display(p_combined)

# =============================================================================
# Summary
# =============================================================================

println("\n" * "=" ^ 70)
println("Summary")
println("=" ^ 70)
println("""

This example demonstrates the new thermodynamic property functions:

1. **thermodynamic_properties(phase, T, y, db)** - Efficient one-call calculation
   Returns NamedTuple: (G, H, S, G_mix, H_mix, S_mix)

2. **Individual functions** (also available):
   - calculate_mixing_gibbs(phase, T, y, db)
   - calculate_mixing_enthalpy(phase, T, y, db)
   - calculate_mixing_entropy(phase, T, y, db)
   - calculate_entropy(phase, T, y, db)
   - calculate_enthalpy(phase, T, y, db)

These functions use automatic differentiation (ForwardDiff) for accurate
temperature derivatives, replacing the numerical differentiation in the
original implementation.

This corresponds to step2.OCM from openCALPHAD (Fortran).
""")
