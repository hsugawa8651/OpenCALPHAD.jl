# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Two-phase equilibrium (map2.OCM equivalent)
#
# Calculate two-phase equilibrium for Cr-Fe system:
# 1. BCC miscibility gap at low temperature
# 2. BCC + LIQUID equilibrium (solidus/liquidus)
#
# This corresponds to map2.OCM from openCALPHAD (Fortran) which calculates
# Cr-Mo miscibility gap and liquidus. We use Cr-Fe here since we have the TDB.
#
# Usage:
#   julia --project examples/312_two_phase_crfe.jl

using OpenCALPHAD
using Plots

println("=" ^ 70)
println("Two-Phase Equilibrium: Cr-Fe System (map2.OCM equivalent)")
println("=" ^ 70)

# Load Cr-Fe thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "crfe.TDB")
db = read_tdb(tdb_path)

bcc = get_phase(db, "BCC_A2")
liquid = get_phase(db, "LIQUID")

println("\nPhases loaded:")
println("  BCC_A2: ", bcc.constituents)
println("  LIQUID: ", liquid.constituents)

# =============================================================================
# Part 1: BCC Miscibility Gap
# =============================================================================

println("\n" * "-" ^ 70)
println("Part 1: BCC Miscibility Gap")
println("-" ^ 70)

# Find miscibility gap at various temperatures
T_gap = Float64[]
x1_gap = Float64[]
x2_gap = Float64[]

println("\nCalculating BCC miscibility gap...")
for T in 500.0:25.0:900.0
    gap = find_miscibility_gap(bcc, T, db)
    if !isnothing(gap) && gap.converged
        push!(T_gap, T)
        push!(x1_gap, gap.x1)  # Cr-rich
        push!(x2_gap, gap.x2)  # Fe-rich
    end
end

if !isempty(T_gap)
    println("Miscibility gap found from T=$(minimum(T_gap))K to T=$(maximum(T_gap))K")
    println(
        "At T=$(T_gap[end])K: x(Cr) = $(round(x1_gap[end], digits=3)) and $(round(x2_gap[end], digits=3))",
    )
else
    println("No miscibility gap found in BCC")
end

# =============================================================================
# Part 2: BCC + LIQUID Equilibrium (Solidus/Liquidus)
# =============================================================================

println("\n" * "-" ^ 70)
println("Part 2: BCC + LIQUID Equilibrium")
println("-" ^ 70)

# Find BCC-LIQUID phase boundaries
T_liq = Float64[]
x_bcc = Float64[]
x_liquid = Float64[]

println("\nCalculating BCC-LIQUID phase boundaries...")
for T in 1600.0:25.0:2200.0
    result = find_common_tangent(bcc, liquid, T, db; n_points = 51)
    if !isnothing(result) && result.converged
        push!(T_liq, T)
        # Get compositions from site fractions
        x_b = result.phase1.site_fractions[1, 1]  # Cr in BCC
        x_l = result.phase2.site_fractions[1, 1]  # Cr in LIQUID
        push!(x_bcc, x_b)
        push!(x_liquid, x_l)
    end
end

if !isempty(T_liq)
    println(
        "BCC-LIQUID equilibrium found from T=$(minimum(T_liq))K to T=$(maximum(T_liq))K",
    )
    println(
        "At T=$(T_liq[1])K: x(Cr)_BCC = $(round(x_bcc[1], digits=3)), x(Cr)_LIQ = $(round(x_liquid[1], digits=3))",
    )
else
    println("No BCC-LIQUID equilibrium found")
end

# =============================================================================
# Create Phase Diagram Plot
# =============================================================================

println("\n" * "-" ^ 70)
println("Generating phase diagram...")
println("-" ^ 70)

p = plot(
    xlabel = "x(Cr)",
    ylabel = "Temperature [K]",
    title = "Cr-Fe Binary Phase Diagram (BCC + LIQUID)",
    xlims = (0, 1),
    ylims = (400, 2200),
    legend = :bottomleft,
    size = (700, 500),
)

# Plot miscibility gap (BCC1 + BCC2)
if !isempty(T_gap)
    plot!(p, x1_gap, T_gap,
        label = "BCC miscibility gap",
        linewidth = 2, color = :blue)
    plot!(p, x2_gap, T_gap,
        label = "",
        linewidth = 2, color = :blue)

    # Shade the two-phase region
    # Fill between the curves
    plot!(p, vcat(x1_gap, reverse(x2_gap)), vcat(T_gap, reverse(T_gap)),
        fill = true, fillalpha = 0.2, fillcolor = :blue,
        linewidth = 0, label = "")
end

# Plot solidus/liquidus (BCC + LIQUID)
if !isempty(T_liq)
    plot!(p, x_bcc, T_liq,
        label = "Solidus (BCC)",
        linewidth = 2, color = :red)
    plot!(p, x_liquid, T_liq,
        label = "Liquidus (LIQUID)",
        linewidth = 2, color = :orange, linestyle = :dash)
end

# Add phase labels
annotate!(p, 0.5, 600, text("BCC₁ + BCC₂", 10, :center))
annotate!(p, 0.5, 1000, text("BCC", 12, :center))
annotate!(p, 0.5, 1900, text("BCC + L", 10, :center))
annotate!(p, 0.5, 2100, text("LIQUID", 12, :center))

# Save
output_path = joinpath(@__DIR__, "312_two_phase_crfe.png")
savefig(p, output_path)
println("\nSaved to: $output_path")

display(p)

# =============================================================================
# Summary
# =============================================================================

println("\n" * "=" ^ 70)
println("Summary")
println("=" ^ 70)
println("""

This example demonstrates two-phase equilibrium calculations:

1. **BCC Miscibility Gap** (low temperature)
   - find_miscibility_gap(bcc, T, db) finds spinodal decomposition
   - Two BCC phases with different Cr content coexist
   - Similar to Ag-Cu FCC miscibility gap

2. **BCC + LIQUID Equilibrium** (high temperature)
   - find_common_tangent(bcc, liquid, T, db) finds solidus/liquidus
   - BCC (solid) and LIQUID phases coexist at melting
   - Uses chemical potential equality condition

This corresponds to map2.OCM from openCALPHAD (Fortran) which calculates
Cr-Mo phase diagram with miscibility gap and liquidus.

Note: The Cr-Fe system also has FCC and SIGMA phases at intermediate
temperatures, which are not shown in this simplified example.
""")
