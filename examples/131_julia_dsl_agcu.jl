# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Julia DSL for Ag-Cu system
#
# Demonstrates how to define thermodynamic models using Julia functions
# instead of TDB files. This enables automatic differentiation through
# the entire calculation.
#
# Usage:
#   julia --project examples/131_julia_dsl_agcu.jl

using OpenCALPHAD

println("=" ^ 60)
println("Julia DSL Example: Ag-Cu FCC Phase")
println("=" ^ 60)

# =============================================================================
# Define phase using Julia DSL (no TDB file needed)
# =============================================================================

# Create FCC phase with single sublattice containing Ag and Cu
fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])

# Set Gibbs energy functions for pure elements (GHSER)
# These are simplified expressions for demonstration
set_G!(fcc, [:AG], T -> -7209.512 + 118.202 * T - 23.8463 * T * log(T))
set_G!(fcc, [:CU], T -> -7770.458 + 130.485 * T - 24.1124 * T * log(T))

# Set interaction parameters (Redlich-Kister)
# L0 and L1 from COST 531 database
set_L!(fcc, [:AG, :CU], 0, T -> 33819.09 - 8.1236 * T)   # L0
set_L!(fcc, [:AG, :CU], 1, T -> -5601.87 + 1.32997 * T)  # L1

println("\nPhase defined:")
println("  Name: $(fcc.name)")
println("  Sublattices: $(length(fcc.sites))")
println("  Constituents: $(fcc.constituents)")

# =============================================================================
# Calculate Gibbs energy at a specific composition
# =============================================================================

println("\n" * "-" ^ 60)
println("Gibbs Energy Calculation")
println("-" ^ 60)

T = 1000.0  # Temperature [K]
x_Cu = 0.3  # Mole fraction of Cu

# Create site fraction matrix
# For single sublattice: y[1,1] = x_Ag, y[1,2] = x_Cu
y = zeros(1, 2)
y[1, 1] = 1.0 - x_Cu  # x_Ag
y[1, 2] = x_Cu        # x_Cu

G = calculate_gibbs_energy(fcc, T, y)

println("\nAt T = $(T) K, x(Cu) = $(x_Cu):")
println("  Gibbs energy: $(round(G, digits=2)) J/mol")

# =============================================================================
# Find miscibility gap
# =============================================================================

println("\n" * "-" ^ 60)
println("Miscibility Gap Calculation")
println("-" ^ 60)

gap = find_miscibility_gap(fcc, T)

if !isnothing(gap)
    println("\nMiscibility gap at T = $(T) K:")
    println(
        "  Ag-rich phase: x(Ag) = $(round(gap.x2, digits=4)), x(Cu) = $(round(1-gap.x2, digits=4))",
    )
    println(
        "  Cu-rich phase: x(Ag) = $(round(gap.x1, digits=4)), x(Cu) = $(round(1-gap.x1, digits=4))",
    )
else
    println("\nNo miscibility gap found at T = $(T) K")
end

# =============================================================================
# Temperature scan
# =============================================================================

println("\n" * "-" ^ 60)
println("Temperature Scan")
println("-" ^ 60)

println("\n  T [K]  |  Ag-rich x(Cu)  |  Cu-rich x(Cu)")
println("-" ^ 50)

for T_scan in [800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0]
    gap_T = find_miscibility_gap(fcc, T_scan)
    if !isnothing(gap_T)
        x_Cu_Ag_rich = 1.0 - gap_T.x2
        x_Cu_Cu_rich = 1.0 - gap_T.x1
        println(
            "  $(Int(T_scan))   |     $(round(x_Cu_Ag_rich, digits=4))      |     $(round(x_Cu_Cu_rich, digits=4))",
        )
    else
        println("  $(Int(T_scan))   |  (no gap)")
    end
end

# =============================================================================
# Advantage: Automatic differentiation works through the calculation
# =============================================================================

println("\n" * "-" ^ 60)
println("Automatic Differentiation")
println("-" ^ 60)

using ForwardDiff

# Define Gibbs energy as a function of composition for AD
function G_of_x(x_Cu_val)
    y_ad = zeros(eltype(x_Cu_val), 1, 2)
    y_ad[1, 1] = 1.0 - x_Cu_val
    y_ad[1, 2] = x_Cu_val
    return calculate_gibbs_energy(fcc, 1000.0, y_ad)
end

# Compute derivative dG/dx using ForwardDiff
x_test = 0.3
dG_dx = ForwardDiff.derivative(G_of_x, x_test)

println("\nAt T = 1000 K, x(Cu) = $(x_test):")
println("  dG/dx(Cu) = $(round(dG_dx, digits=2)) J/mol")

println("\n" * "=" ^ 60)
println("Example completed successfully!")
println("=" ^ 60)
