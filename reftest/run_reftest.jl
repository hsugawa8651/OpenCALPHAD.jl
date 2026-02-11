#!/usr/bin/env julia
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Reference test runner
#
# Reference test runner for OpenCALPHAD.jl
# Compares Julia implementation results with original openCALPHAD
#

using JSON
using Test

# Add parent directory to load path
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using OpenCALPHAD

const REFTEST_DIR = @__DIR__
const TDB_DIR = joinpath(REFTEST_DIR, "tdb")
const REF_DIR = joinpath(REFTEST_DIR, "reference")

"""
    load_reference(name::String) -> Dict

Load reference data from JSON file.
"""
function load_reference(name::String)
    path = joinpath(REF_DIR, name * ".json")
    if !isfile(path)
        error("Reference file not found: $path")
    end
    return JSON.parsefile(path)
end

"""
    ReftestResult

Result of a reference test comparison.
"""
struct ReftestResult
    name::String
    passed::Bool
    errors::Vector{String}
    computed::Dict{String, Any}
    reference::Dict{String, Any}
end

"""
    compare_gibbs_energy(computed, reference; atol=10.0) -> Vector{String}

Compare total Gibbs energy. Default tolerance: ±10 J/mol.
"""
function compare_gibbs_energy(computed::Float64, reference::Float64; atol::Float64 = 10.0)
    errors = String[]
    if !isapprox(computed, reference, atol = atol)
        push!(
            errors,
            "G mismatch: computed=$(round(computed, digits=1)), reference=$reference, diff=$(round(computed-reference, digits=1)) J/mol",
        )
    end
    return errors
end

"""
    compare_phase_amounts(computed, reference; atol=0.01) -> Vector{String}

Compare phase amounts. Default tolerance: ±0.01 mol.
"""
function compare_phase_amounts(
    computed_phases::Vector,
    ref_phases::Vector;
    atol::Float64 = 0.01,
)
    errors = String[]

    for ref_phase in ref_phases
        ref_name = ref_phase["name"]
        ref_amount = ref_phase["amount"]

        # Find matching computed phase (by similar composition)
        found = false
        for (i, comp_phase) in enumerate(computed_phases)
            comp_amount = comp_phase["amount"]
            if isapprox(comp_amount, ref_amount, atol = atol)
                found = true
                break
            end
        end

        if !found
            push!(errors, "Phase amount mismatch for $ref_name: expected $ref_amount")
        end
    end

    return errors
end

"""
    compare_compositions(computed, reference; atol=0.01) -> Vector{String}

Compare phase compositions. Default tolerance: ±0.01 mole fraction.
"""
function compare_compositions(
    computed_phases::Vector,
    ref_phases::Vector;
    atol::Float64 = 0.01,
)
    errors = String[]

    # Match phases by AG composition (more reliable than Dict iteration order)
    for ref_phase in ref_phases
        ref_comp = ref_phase["composition"]
        ref_x_ag = get(ref_comp, "AG", 0.0)

        matched = false
        for comp_phase in computed_phases
            comp_comp = comp_phase["composition"]
            comp_x_ag = get(comp_comp, "AG", 0.0)

            # Match by AG composition
            if isapprox(comp_x_ag, ref_x_ag, atol = atol)
                matched = true
                # Check all components
                for (elem, ref_val) in ref_comp
                    comp_val = get(comp_comp, uppercase(string(elem)), NaN)
                    if !isapprox(comp_val, ref_val, atol = atol)
                        push!(
                            errors,
                            "Composition mismatch: $elem in $(ref_phase["name"]): computed=$(round(comp_val, digits=4)), reference=$ref_val",
                        )
                    end
                end
                break
            end
        end

        if !matched
            push!(
                errors,
                "No matching phase found for $(ref_phase["name"]) with x(AG)=$(ref_x_ag)",
            )
        end
    end

    return errors
end

"""
    run_agcu_equil_1000K() -> ReftestResult

Run Ag-Cu equilibrium test at T=1000K, x(Cu)=0.3.
"""
function run_agcu_equil_1000K()
    name = "agcu_equil_1000K"
    ref = load_reference(name)

    # Load database
    tdb_path = joinpath(TDB_DIR, "agcu.TDB")
    db = read_tdb(tdb_path)

    # Get FCC phase
    fcc = get_phase(db, "FCC_A1")

    # Conditions
    T = Float64(ref["conditions"]["T"])
    x_Ag = 1.0 - Float64(ref["conditions"]["x_Cu"])  # x_Ag = 0.7

    # Find miscibility gap
    gap = find_miscibility_gap(fcc, T, db, 0.035, 0.90)

    if isnothing(gap)
        return ReftestResult(name, false, ["Miscibility gap not found"], Dict(), ref)
    end

    # Calculate equilibrium at x_Ag = 0.7
    result = calculate_miscibility_equilibrium(fcc, T, db, x_Ag)

    if isnothing(result)
        return ReftestResult(name, false, ["Equilibrium calculation failed"], Dict(), ref)
    end

    # Extract computed values
    G_total =
        result.phase1_fraction * result.phase1.gibbs_energy +
        result.phase2_fraction * result.phase2.gibbs_energy

    computed_phases = [
        Dict(
            "name" => "FCC_A1#1",
            "amount" => result.phase2_fraction,  # Ag-rich
            "composition" => Dict(
                "AG" => result.phase2.site_fractions[1, 1],
                "CU" => result.phase2.site_fractions[1, 2],
            ),
        ),
        Dict(
            "name" => "FCC_A1#2",
            "amount" => result.phase1_fraction,  # Cu-rich
            "composition" => Dict(
                "AG" => result.phase1.site_fractions[1, 1],
                "CU" => result.phase1.site_fractions[1, 2],
            ),
        ),
    ]

    computed = Dict{String, Any}(
        "G_total" => G_total,
        "phases" => computed_phases,
    )

    # Compare results
    errors = String[]

    # Compare G
    ref_G = Float64(ref["results"]["G_total"])
    append!(errors, compare_gibbs_energy(G_total, ref_G))

    # Compare phase amounts
    ref_phases = ref["results"]["phases"]
    append!(errors, compare_phase_amounts(computed_phases, ref_phases))

    # Compare compositions
    append!(errors, compare_compositions(computed_phases, ref_phases))

    passed = isempty(errors)
    return ReftestResult(name, passed, errors, computed, ref)
end

"""
    run_all_reftests() -> Vector{ReftestResult}

Run all reference tests.
"""
function run_all_reftests()
    results = ReftestResult[]

    # Add all test functions here
    push!(results, run_agcu_equil_1000K())

    return results
end

"""
    print_results(results::Vector{ReftestResult})

Print test results summary.
"""
function print_results(results::Vector{ReftestResult})
    println("\n" * "="^60)
    println("Reference Test Results")
    println("="^60)

    passed = 0
    failed = 0

    for r in results
        status = r.passed ? "PASS" : "FAIL"
        println("\n[$status] $(r.name)")

        if r.passed
            passed += 1
            if haskey(r.computed, "G_total")
                println("  G = $(round(r.computed["G_total"], digits=1)) J/mol")
            end
        else
            failed += 1
            for err in r.errors
                println("  ERROR: $err")
            end
        end
    end

    println("\n" * "-"^60)
    println("Summary: $passed passed, $failed failed")
    println("="^60)

    return failed == 0
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("OpenCALPHAD.jl Reference Tests")
    println("Comparing with openCALPHAD v6.100")
    println()

    results = run_all_reftests()
    success = print_results(results)

    exit(success ? 0 : 1)
end

# Export for use in test suite
export run_all_reftests, run_agcu_equil_1000K, ReftestResult
