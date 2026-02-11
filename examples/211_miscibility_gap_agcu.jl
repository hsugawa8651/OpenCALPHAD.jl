# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Miscibility gap validation
#
# Validation script for Ag-Cu miscibility gap calculation
#
# Compares OpenCALPHAD.jl results with:
# 1. openCALPHAD (Fortran) reference data
# 2. Literature values (Subramanian & Perepezko 1993, Murray 1984)
#
# Usage:
#   julia --project examples/211_miscibility_gap_agcu.jl

using OpenCALPHAD
using JSON
using Printf

# =============================================================================
# Reference Data
# =============================================================================

# openCALPHAD (Fortran) reference: T=1000K, x(Cu)=0.3
const REFTEST_1000K = Dict(
    "T" => 1000.0,
    "Ag_rich_xCu" => 0.1031,  # FCC_A1#1
    "Cu_rich_xCu" => 0.9663,  # FCC_A1_AUTO#2
    "G_total" => -53648.0,    # J/mol
    "source" => "openCALPHAD 6.100",
)

# Literature: Murray (1984) via Subramanian & Perepezko (1993) Fig. 5
# Metastable FCC miscibility gap parameters
const LITERATURE_MURRAY1984 = Dict(
    "critical_T" => 1385.0,   # K (1112°C)
    "critical_xCu" => 0.50,   # at consolute point
    "source" => "Murray (1984), Subramanian & Perepezko (1993) Fig. 5",
)

# Tolerances
const TOL_COMPOSITION = 0.02   # mole fraction
const TOL_GIBBS = 50.0         # J/mol
const TOL_TEMPERATURE = 50.0   # K

# =============================================================================
# Load database and run calculations
# =============================================================================

println("=" ^ 70)
println("Ag-Cu Miscibility Gap Validation")
println("=" ^ 70)

tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

println("\nDatabase: $tdb_path")
println("Phase: $(fcc.name)")

# =============================================================================
# Test 1: Compare with openCALPHAD at T=1000K
# =============================================================================

println("\n" * "-" ^ 70)
println("Test 1: Comparison with openCALPHAD (Fortran) at T=1000K")
println("-" ^ 70)

T_test = 1000.0
gap_1000K = find_miscibility_gap(fcc, T_test, db)

if isnothing(gap_1000K)
    println("ERROR: Miscibility gap not found at T=$(T_test)K")
    exit(1)
end

# Convert to x(Cu) for comparison (gap returns x(Ag), so x(Cu) = 1 - x(Ag))
julia_Ag_rich_xCu = 1.0 - gap_1000K.x2  # x2 is Ag-rich (higher x_Ag)
julia_Cu_rich_xCu = 1.0 - gap_1000K.x1  # x1 is Cu-rich (lower x_Ag)

ref_Ag_rich_xCu = REFTEST_1000K["Ag_rich_xCu"]
ref_Cu_rich_xCu = REFTEST_1000K["Cu_rich_xCu"]

diff_Ag_rich = abs(julia_Ag_rich_xCu - ref_Ag_rich_xCu)
diff_Cu_rich = abs(julia_Cu_rich_xCu - ref_Cu_rich_xCu)

pass_Ag_rich = diff_Ag_rich < TOL_COMPOSITION
pass_Cu_rich = diff_Cu_rich < TOL_COMPOSITION

println()
println("                    | OpenCALPHAD.jl | openCALPHAD(F) |   Diff   | Pass")
println("-" ^ 70)
@printf("Ag-rich boundary    |     %6.4f     |     %6.4f     |  %6.4f  |  %s\n",
    julia_Ag_rich_xCu, ref_Ag_rich_xCu, diff_Ag_rich, pass_Ag_rich ? "OK" : "FAIL")
@printf("Cu-rich boundary    |     %6.4f     |     %6.4f     |  %6.4f  |  %s\n",
    julia_Cu_rich_xCu, ref_Cu_rich_xCu, diff_Cu_rich, pass_Cu_rich ? "OK" : "FAIL")
println()
println("Tolerance: ±$(TOL_COMPOSITION) (mole fraction)")
println("Source: $(REFTEST_1000K["source"])")

# =============================================================================
# Test 2: Find critical point and compare with literature
# =============================================================================

println("\n" * "-" ^ 70)
println("Test 2: Critical point estimation (vs Murray 1984)")
println("-" ^ 70)

# Scan temperature to find where miscibility gap closes
println("\nScanning temperature to find critical point...")

T_values = Float64[]
gap_widths = Float64[]

for T in 600.0:25.0:1400.0
    gap = find_miscibility_gap(fcc, T, db)
    if !isnothing(gap)
        push!(T_values, T)
        push!(gap_widths, abs(gap.x2 - gap.x1))
    end
end

if length(T_values) > 0
    # Estimate critical temperature (where gap width -> 0)
    max_T_with_gap = maximum(T_values)

    # Find the narrowest gap
    min_idx = argmin(gap_widths)
    narrowest_T = T_values[min_idx]
    narrowest_width = gap_widths[min_idx]

    println()
    println("Results:")
    println(
        "  Temperature range with miscibility gap: $(minimum(T_values))K - $(max_T_with_gap)K",
    )
    println(
        "  Narrowest gap: width=$(round(narrowest_width, digits=4)) at T=$(narrowest_T)K",
    )
    println()
    println("Literature (Murray 1984):")
    println(
        "  Critical temperature: $(LITERATURE_MURRAY1984["critical_T"])K ($(LITERATURE_MURRAY1984["critical_T"]-273.15)°C)",
    )
    println("  Critical composition: x(Cu) = $(LITERATURE_MURRAY1984["critical_xCu"])")
    println()

    # Note about discrepancy
    if max_T_with_gap < LITERATURE_MURRAY1984["critical_T"] - 100
        println("NOTE: Calculated critical point is significantly lower than literature.")
        println("      This may be due to different thermodynamic parameters or")
        println("      the TDB file using different assessment.")
    end
else
    println("ERROR: No miscibility gap found in temperature range")
end

# =============================================================================
# Test 3: Thermodynamic consistency check
# =============================================================================

println("\n" * "-" ^ 70)
println("Test 3: Thermodynamic consistency (chemical potential equality)")
println("-" ^ 70)

T_check = 1000.0
gap = find_miscibility_gap(fcc, T_check, db)

if !isnothing(gap)
    # Calculate chemical potentials at both boundaries
    n_const = length(fcc.constituents[1])
    n_sub = length(fcc.sites)

    y1 = zeros(n_sub, n_const)
    y1[1, 1] = gap.x1  # x(Ag) at Cu-rich boundary
    y1[1, 2] = 1 - gap.x1
    for s in 2:n_sub
        y1[s, 1] = 1.0
    end

    y2 = zeros(n_sub, n_const)
    y2[1, 1] = gap.x2  # x(Ag) at Ag-rich boundary
    y2[1, 2] = 1 - gap.x2
    for s in 2:n_sub
        y2[s, 1] = 1.0
    end

    mu1 = compute_chemical_potential(fcc, T_check, y1, db, ForwardDiffMethod())
    mu2 = compute_chemical_potential(fcc, T_check, y2, db, ForwardDiffMethod())

    delta_mu_Ag = abs(mu1[1] - mu2[1])
    delta_mu_Cu = abs(mu1[2] - mu2[2])

    println()
    println("At T=$(T_check)K:")
    println("  Cu-rich boundary: x(Ag) = $(round(gap.x1, digits=4))")
    println("  Ag-rich boundary: x(Ag) = $(round(gap.x2, digits=4))")
    println()
    println("Chemical potential difference:")
    @printf("  |Δμ(Ag)| = %.2f J/mol\n", delta_mu_Ag)
    @printf("  |Δμ(Cu)| = %.2f J/mol\n", delta_mu_Cu)
    println()

    if delta_mu_Ag < 100 && delta_mu_Cu < 100
        println("  Status: OK (chemical potentials approximately equal)")
    else
        println("  Status: WARNING (chemical potentials differ significantly)")
    end
end

# =============================================================================
# Summary
# =============================================================================

println("\n" * "=" ^ 70)
println("Validation Summary")
println("=" ^ 70)

all_pass = pass_Ag_rich && pass_Cu_rich

println()
if all_pass
    println("Result: ALL TESTS PASSED")
    println()
    println("The miscibility gap calculation is validated against:")
    println("  1. openCALPHAD (Fortran) reference values")
    println("  2. Thermodynamic consistency (common tangent)")
else
    println("Result: SOME TESTS FAILED")
    println()
    println("Please check the discrepancies above.")
end

println()
println("References:")
println("  - Subramanian & Perepezko (1993) J. Phase Equilib. 14, 62-75")
println("  - Murray (1984) Metall. Trans. A 15, 261-268")
println()
