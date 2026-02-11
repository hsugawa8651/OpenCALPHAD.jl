# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Validated phase diagram with Julia DSL
#
# Ag-Cu Phase Diagram with Validation (Julia DSL)
#
# Combined figure showing:
# - OpenCALPHAD.jl calculated phase boundaries using Julia DSL (lines)
# - openCALPHAD (Fortran) reference data (markers)
#
# This demonstrates that the Julia DSL produces identical results to TDB.
#
# Usage:
#   julia --project examples/911_phase_diagram_validated.jl

using OpenCALPHAD
using Plots
using Printf

println("=" ^ 70)
println("Ag-Cu Phase Diagram with Validation (Julia DSL)")
println("=" ^ 70)

# =============================================================================
# Reference data from openCALPHAD (Fortran) v6.100
# =============================================================================

const FORTRAN_REF = Dict(
    # T => (x_Ag_rich, x_Cu_rich) where x is x(Cu)
    800.0 => (0.0380, 0.9910),
    900.0 => (0.0651, 0.9817),
    1000.0 => (0.1031, 0.9663),
    1100.0 => (0.1545, 0.9424),
    1200.0 => (0.2240, 0.9052),
)

# =============================================================================
# Define Ag-Cu FCC phase using Julia DSL
# =============================================================================

println("\nDefining phase using Julia DSL...")

# GHSER functions for pure elements (from COST 531 database)
function ghser_ag(T)
    return -7209.512 + 118.200733*T - 23.8463314*T*log(T) -
           0.001790585*T^2 - 3.98587e-7*T^3 - 12011/T
end

function ghser_cu(T)
    return -7770.458 + 130.485403*T - 24.112392*T*log(T) -
           0.00265684*T^2 + 1.29223e-7*T^3 + 52478/T
end

# Create FCC phase with Julia DSL
fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])

set_G!(fcc, [:AG], ghser_ag)              # G of pure Ag
set_G!(fcc, [:CU], ghser_cu)              # G of pure Cu
set_L!(fcc, [:AG, :CU], 0, T -> 33819.09 - 8.1236*T)   # L₀, 0th order
set_L!(fcc, [:AG, :CU], 1, T -> -5601.87 + 1.32997*T)  # L₁, 1st order

println("  Phase defined: $(fcc.name)")
println("  Julia parameters: $(length(fcc.julia_parameters))")

# =============================================================================
# Calculate phase boundaries
# =============================================================================

println("\nCalculating phase boundaries (800-1300 K)...")

T_range = 800.0:25.0:1300.0
julia_T = Float64[]
julia_x1 = Float64[]  # Ag-rich boundary (lower x_Cu)
julia_x2 = Float64[]  # Cu-rich boundary (higher x_Cu)

for T in T_range
    gap = find_miscibility_gap(fcc, T)  # DSL version (no db argument)
    if !isnothing(gap)
        push!(julia_T, T)
        # gap.x1, gap.x2 are x(Ag), convert to x(Cu) = 1 - x(Ag)
        push!(julia_x1, 1.0 - gap.x2)  # Ag-rich: higher x(Ag) → lower x(Cu)
        push!(julia_x2, 1.0 - gap.x1)  # Cu-rich: lower x(Ag) → higher x(Cu)
    end
end

println("  Calculated $(length(julia_T)) temperature points")

# =============================================================================
# Create combined figure
# =============================================================================

println("\nCreating figure...")

# Extract reference data
ref_T = sort(collect(keys(FORTRAN_REF)))
ref_x1 = [FORTRAN_REF[T][1] for T in ref_T]
ref_x2 = [FORTRAN_REF[T][2] for T in ref_T]

p = plot(
    xlabel = "x(Cu)",
    ylabel = "Temperature [K]",
    title = "Ag-Cu FCC Miscibility Gap",
    legend = :top,
    size = (650, 500),
    xlims = (0, 1),
    ylims = (750, 1350),
    grid = true,
    framestyle = :box,
)

# Fill two-phase region
x_fill = vcat(julia_x1, reverse(julia_x2))
T_fill = vcat(julia_T, reverse(julia_T))
plot!(p, x_fill, T_fill,
    fillrange = minimum(julia_T),
    fillalpha = 0.15,
    fillcolor = :blue,
    linewidth = 0,
    label = "",
)

# OpenCALPHAD.jl calculated lines (Julia DSL)
plot!(p, julia_x1, julia_T,
    label = "OpenCALPHAD.jl (Julia DSL)",
    linewidth = 2.5,
    color = :blue,
)
plot!(p, julia_x2, julia_T,
    label = "",
    linewidth = 2.5,
    color = :blue,
)

# openCALPHAD (Fortran) reference points
scatter!(p, ref_x1, ref_T,
    label = "openCALPHAD (Fortran)",
    markersize = 8,
    color = :red,
    markershape = :circle,
    markerstrokewidth = 2,
    markerstrokecolor = :red,
)
scatter!(p, ref_x2, ref_T,
    label = "",
    markersize = 8,
    color = :red,
    markershape = :circle,
    markerstrokewidth = 2,
    markerstrokecolor = :red,
)

# Annotations
annotate!(p, 0.5, 900, text("Two-phase\nregion", 10, :center))

# Save
output_path = joinpath(@__DIR__, "911_phase_diagram_validated.png")
savefig(p, output_path)
println("Saved to: $output_path")

# =============================================================================
# Print numerical comparison
# =============================================================================

println("\n" * "-" ^ 60)
println("Validation: OpenCALPHAD.jl (Julia DSL) vs openCALPHAD (Fortran)")
println("-" ^ 60)
println()
println("T [K]  | Ag-rich x(Cu)      | Cu-rich x(Cu)      |")
println("       | Fortran  Julia     | Fortran  Julia     |")
println("-" ^ 60)

max_diff = 0.0
for T in ref_T
    idx = argmin(abs.(julia_T .- T))
    if abs(julia_T[idx] - T) < 1.0
        j_x1 = julia_x1[idx]
        j_x2 = julia_x2[idx]
        r_x1 = FORTRAN_REF[T][1]
        r_x2 = FORTRAN_REF[T][2]
        diff1 = abs(j_x1 - r_x1)
        diff2 = abs(j_x2 - r_x2)
        global max_diff = max(max_diff, diff1, diff2)
        @printf("%6.0f | %7.4f  %7.4f   | %7.4f  %7.4f   |\n",
            T, r_x1, j_x1, r_x2, j_x2)
    end
end

println("-" ^ 60)
@printf("Maximum difference: %.4f (%.2f%%)\n", max_diff, max_diff * 100)
println()

display(p)
