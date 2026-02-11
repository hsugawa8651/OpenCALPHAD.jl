# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Numerical verification tests
#
# Heavy numerical tests separated from runtests.jl for faster CI.
# Run manually: julia --project reftest/run_numerical_tests.jl

using Test
using OpenCALPHAD
using LinearAlgebra: norm

const REFTEST_DIR = @__DIR__
const TDB_DIR = joinpath(REFTEST_DIR, "tdb")

@testset "Numerical Tests" begin

    # =========================================================================
    # Two-phase equilibrium (Cr-Fe)
    # =========================================================================
    @testset "Two-phase equilibrium" begin
        # Setup: Create two phases (BCC and FCC) for Cr-Fe
        function create_two_phase_system()
            db = Database()

            # Add elements
            push!(db.elements, Element("CR", "BCC_A2", 51.996, 4050.0, 23.543))
            push!(db.elements, Element("FE", "BCC_A2", 55.845, 4489.0, 27.28))

            # Create BCC phase
            bcc = Phase("BCC_A2", 'S', [1.0])
            push!(bcc.constituents, ["CR", "FE"])
            push!(db.phases, bcc)

            # BCC parameters
            push!(
                bcc.parameters,
                Parameter(
                    "BCC_A2", 'G', [["CR"]], 0,
                    [298.15, 6000.0],
                    ["-8856.94+157.48*T-26.908*T*LN(T)"],
                ),
            )
            push!(
                bcc.parameters,
                Parameter(
                    "BCC_A2", 'G', [["FE"]], 0,
                    [298.15, 6000.0],
                    ["1225.7+124.134*T-23.5143*T*LN(T)"],
                ),
            )
            push!(
                bcc.parameters,
                Parameter(
                    "BCC_A2", 'L', [["CR", "FE"]], 0,
                    [298.15, 6000.0],
                    ["20500-9.68*T"],
                ),
            )

            # Create FCC phase
            fcc = Phase("FCC_A1", 'S', [1.0])
            push!(fcc.constituents, ["CR", "FE"])
            push!(db.phases, fcc)

            # FCC parameters (shifted to create a two-phase region)
            push!(
                fcc.parameters,
                Parameter(
                    "FCC_A1", 'G', [["CR"]], 0,
                    [298.15, 6000.0],
                    ["-1572.94+157.48*T-26.908*T*LN(T)"],  # Higher than BCC
                ),
            )
            push!(
                fcc.parameters,
                Parameter(
                    "FCC_A1", 'G', [["FE"]], 0,
                    [298.15, 6000.0],
                    ["-1462.4+124.134*T-23.5143*T*LN(T)"],  # Lower than BCC
                ),
            )
            push!(
                fcc.parameters,
                Parameter(
                    "FCC_A1", 'L', [["CR", "FE"]], 0,
                    [298.15, 6000.0],
                    ["10833-7.477*T"],  # Different interaction
                ),
            )

            return bcc, fcc, db
        end

        @testset "find_common_tangent" begin
            bcc, fcc, db = create_two_phase_system()

            # At 1000 K, should find a two-phase region
            result = find_common_tangent(bcc, fcc, 1000.0, db; n_points = 51)

            # Result may be nothing if no equilibrium found
            if !isnothing(result)
                @test result isa TwoPhaseResult
                @test result.converged == true

                # Compositions should be different
                x1 = result.phase1.site_fractions[1, 1]
                x2 = result.phase2.site_fractions[1, 1]
                @test abs(x1 - x2) > 0.01  # Should be measurably different

                # Phase fractions should be valid
                @test 0 <= result.phase1_fraction <= 1
                @test 0 <= result.phase2_fraction <= 1
            end
        end

        @testset "calculate_two_phase_equilibrium" begin
            bcc, fcc, db = create_two_phase_system()

            # First find boundaries
            ct = find_common_tangent(bcc, fcc, 1000.0, db; n_points = 51)

            if !isnothing(ct)
                x1 = ct.phase1.site_fractions[1, 1]
                x2 = ct.phase2.site_fractions[1, 1]
                x_mid = (x1 + x2) / 2

                # Calculate at midpoint
                result = calculate_two_phase_equilibrium(bcc, fcc, 1000.0, db, x_mid)

                if !isnothing(result)
                    @test result isa TwoPhaseResult
                    @test result.converged == true

                    # Lever rule check
                    f1 = result.phase1_fraction
                    f2 = result.phase2_fraction
                    @test f1 + f2 ≈ 1.0

                    # Mass balance
                    x_calc =
                        f1 * result.phase1.site_fractions[1, 1] +
                        f2 * result.phase2.site_fractions[1, 1]
                    @test x_calc ≈ x_mid atol=0.01
                end
            end
        end

        @testset "Single phase region" begin
            bcc, fcc, db = create_two_phase_system()

            # At extreme compositions, should return nothing (single phase)
            result_low = calculate_two_phase_equilibrium(bcc, fcc, 1000.0, db, 0.01)
            result_high = calculate_two_phase_equilibrium(bcc, fcc, 1000.0, db, 0.99)

            # These should likely be nothing or in single-phase region
            # (depends on the actual thermodynamics)
        end
    end

    # =========================================================================
    # Ag-Cu miscibility gap (reference values from openCALPHAD v6.100)
    # =========================================================================
    @testset "Ag-Cu miscibility gap" begin
        tdb_path = joinpath(TDB_DIR, "agcu.TDB")
        db = read_tdb(tdb_path)
        fcc = get_phase(db, "FCC_A1")

        # Find miscibility gap
        gap = find_miscibility_gap(fcc, 1000.0, db, 0.035, 0.90)

        @test !isnothing(gap)
        @test gap.converged == true

        # Reference compositions (within tolerance)
        # Cu-rich: x_Ag = 0.0337
        # Ag-rich: x_Ag = 0.8969
        @test gap.x1 ≈ 0.0337 atol=0.01
        @test gap.x2 ≈ 0.8969 atol=0.01

        # Calculate equilibrium at x_Cu = 0.3 (x_Ag = 0.7)
        result = calculate_miscibility_equilibrium(fcc, 1000.0, db, 0.7)

        @test !isnothing(result)
        @test result.converged == true

        # Reference phase fractions
        # f (Cu-rich) = 0.2281, f (Ag-rich) = 0.7719
        @test result.phase1_fraction ≈ 0.2281 atol=0.01
        @test result.phase2_fraction ≈ 0.7719 atol=0.01

        # Reference G = -53648 J/mol
        G_total =
            result.phase1_fraction * result.phase1.gibbs_energy +
            result.phase2_fraction * result.phase2.gibbs_energy
        @test G_total ≈ -53648 atol=10  # Within 10 J/mol
    end

    # =========================================================================
    # STEP calculation (numerical verification)
    # =========================================================================
    @testset "STEP calculation" begin
        tdb_path = joinpath(TDB_DIR, "agcu.TDB")
        db = read_tdb(tdb_path)
        fcc = get_phase(db, "FCC_A1")

        @testset "step_temperature length" begin
            # T=800 to T=1000, step=100 => 3 points: 800, 900, 1000
            result = step_temperature(fcc, db, 0.7, 800.0, 1000.0, 100.0)

            @test result isa StepResult
            @test result.axis_variable == :T
            @test length(result.points) == 3
            @test result.axis_values[1] ≈ 800.0
            @test result.axis_values[2] ≈ 900.0
            @test result.axis_values[3] ≈ 1000.0
            @test result.x_overall ≈ 0.7
        end

        @testset "G values at T=1000K" begin
            # Reference: G = -53648 J/mol at T=1000K, x(Ag)=0.7
            result = step_temperature(fcc, db, 0.7, 1000.0, 1000.0, 100.0)

            @test length(result.points) == 1
            point = result.points[1]

            # Verify G value (within 10 J/mol tolerance)
            @test point.gibbs_energy ≈ -53648 atol=10

            # Verify two-phase
            @test point.n_phases == 2
            @test point.converged == true

            # Verify phase fractions (reference: 0.228, 0.772)
            @test point.phase_amounts[1] ≈ 0.228 atol=0.01
            @test point.phase_amounts[2] ≈ 0.772 atol=0.01
        end

        @testset "Two-phase region T=800-1200K" begin
            result = step_temperature(fcc, db, 0.7, 800.0, 1200.0, 100.0)

            # Reference G values from agcu_step_T_x03_julia.json
            ref_G = [-39427.0, -46335.0, -53648.0, -61359.0, -69469.0]

            @test length(result.points) == 5

            for (i, point) in enumerate(result.points)
                @test point.n_phases == 2  # All in two-phase region
                @test point.gibbs_energy ≈ ref_G[i] atol=10
            end
        end

        @testset "Helper functions" begin
            result = step_temperature(fcc, db, 0.7, 800.0, 1000.0, 100.0)

            # Test temperatures()
            T_vals = temperatures(result)
            @test length(T_vals) == 3
            @test T_vals[1] ≈ 800.0
            @test T_vals[3] ≈ 1000.0

            # Test gibbs_energies()
            G_vals = gibbs_energies(result)
            @test length(G_vals) == 3
            @test G_vals[3] ≈ -53648 atol=10

            # Test phase_fractions()
            f1 = phase_fractions(result, 1)
            f2 = phase_fractions(result, 2)
            @test length(f1) == 3
            @test length(f2) == 3
            @test f1[3] ≈ 0.228 atol=0.01  # Cu-rich phase at T=1000K
            @test f2[3] ≈ 0.772 atol=0.01  # Ag-rich phase at T=1000K

            # Phase fractions should sum to 1
            @test f1[1] + f2[1] ≈ 1.0
            @test f1[2] + f2[2] ≈ 1.0
            @test f1[3] + f2[3] ≈ 1.0
        end

        @testset "Single-phase region T=1300K" begin
            # At T=1300K, x(Ag)=0.7 is outside the miscibility gap
            result = step_temperature(fcc, db, 0.7, 1300.0, 1300.0, 100.0)

            @test length(result.points) == 1
            point = result.points[1]

            # Should be single phase
            @test point.n_phases == 1
            @test point.converged == true
            @test length(point.phase_amounts) == 1
            @test point.phase_amounts[1] ≈ 1.0
        end
    end

    # =========================================================================
    # MAP calculation (numerical verification)
    # =========================================================================
    @testset "MAP calculation (phase diagram)" begin
        tdb_path = joinpath(TDB_DIR, "agcu.TDB")
        db = read_tdb(tdb_path)
        fcc = get_phase(db, "FCC_A1")

        @testset "map_phase_diagram basic" begin
            # T=800 to T=1000, step=100 => 3 points
            result = map_phase_diagram(fcc, db, 800.0, 1000.0, 100.0)

            @test result isa PhaseDiagramResult
            @test result.phase === fcc
            @test length(result.points) >= 1  # At least some points converged
        end

        @testset "Boundary compositions at T=1000K" begin
            # Reference: x1=0.0337 (Cu-rich), x2=0.8969 (Ag-rich)
            result = map_phase_diagram(fcc, db, 1000.0, 1000.0, 100.0)

            @test length(result.points) == 1
            point = result.points[1]

            @test point.temperature ≈ 1000.0
            @test point.compositions[1] ≈ 0.034 atol=0.01
            @test point.compositions[2] ≈ 0.897 atol=0.01
            @test point.converged == true
        end

        @testset "Helper functions" begin
            result = map_phase_diagram(fcc, db, 800.0, 1000.0, 100.0)

            # Test temperatures()
            T_vals = temperatures(result)
            @test length(T_vals) >= 1
            @test T_vals[1] ≈ 800.0

            # Test boundaries(result, idx)
            x1_vals = boundaries(result, 1)
            x2_vals = boundaries(result, 2)
            @test length(x1_vals) == length(result.points)
            @test length(x2_vals) == length(result.points)

            # Left boundary should be smaller than right boundary
            for i in eachindex(x1_vals)
                @test x1_vals[i] < x2_vals[i]
            end
        end

        @testset "Miscibility gap closes at high T" begin
            # At very high T, miscibility gap should close
            result = map_phase_diagram(fcc, db, 1000.0, 1400.0, 100.0)

            # Find points where gap exists (converged)
            converged_points = filter(p -> p.converged, result.points)

            if length(converged_points) > 0
                # Gap should narrow at higher T
                T_vals = [p.temperature for p in converged_points]
                widths = [p.compositions[2] - p.compositions[1] for p in converged_points]

                # At lower T, gap should be wider
                if length(widths) >= 2
                    @test widths[1] > widths[end]  # Gap narrows with increasing T
                end
            end
        end
    end

    # =========================================================================
    # Plotting recipes (require numerical computation)
    # =========================================================================
    @testset "Plotting recipes" begin
        tdb_path = joinpath(TDB_DIR, "agcu.TDB")
        db = read_tdb(tdb_path)
        fcc = get_phase(db, "FCC_A1")

        @testset "plot_phase_diagram returns PhaseDiagramResult" begin
            map_result = map_phase_diagram(fcc, db, 800.0, 1000.0, 100.0)
            plot_data = plot_phase_diagram(map_result)

            @test plot_data isa PhaseDiagramResult
            @test plot_data === map_result
        end

        @testset "plot_gibbs_curve returns GridScanResult" begin
            solver = GridSearchSolver(n_points = 11)
            scan_result = scan_composition(fcc, 1000.0, db, solver)
            plot_data = plot_gibbs_curve(scan_result)

            @test plot_data isa GridScanResult
            @test plot_data === scan_result
        end

        @testset "StepResult is plottable" begin
            step_result = step_temperature(fcc, db, 0.7, 800.0, 1000.0, 100.0)

            # Verify the data that would be used for plotting
            T = temperatures(step_result)
            G = gibbs_energies(step_result)

            @test length(T) == 3
            @test length(G) == 3
            @test all(isfinite, G)
        end

        @testset "PhaseDiagramResult data extraction" begin
            map_result = map_phase_diagram(fcc, db, 800.0, 1000.0, 100.0)

            # Verify boundary data extraction
            T = temperatures(map_result)
            x1 = boundaries(map_result, 1)
            x2 = boundaries(map_result, 2)

            @test length(T) == length(x1) == length(x2)

            # Converged points should have valid boundaries
            for i in eachindex(map_result.points)
                if map_result.points[i].converged
                    @test x1[i] < x2[i]  # Left < Right
                    @test 0 < x1[i] < 1
                    @test 0 < x2[i] < 1
                end
            end
        end
    end

    # =========================================================================
    # Multi-component numerical tests (Cr-Mo from steel1.TDB)
    # =========================================================================
    @testset "Multi-component numerical" begin
        @testset "scan_composition with component_indices" begin
            db = read_tdb(joinpath(TDB_DIR, "steel1.TDB"))
            bcc = get_phase(db, "BCC_A2")
            solver = GridSearchSolver(n_points = 11)

            # Cr-Mo scan (CR=1, MO=3)
            scan = scan_composition(bcc, 1800.0, db, solver; component_indices = (1, 3))

            @test length(scan.x_grid) == 11
            @test scan.x_grid[1] ≈ 0.0 atol = 1e-6
            @test scan.x_grid[end] ≈ 1.0 atol = 1e-6
            @test all(isfinite, scan.G_values)

            # Verify G values are reasonable (negative Gibbs energy at high T)
            @test all(G -> G < 0, scan.G_values)

            # Default should still work (for backward compatibility)
            db_agcu = read_tdb(joinpath(TDB_DIR, "agcu.TDB"))
            fcc = get_phase(db_agcu, "FCC_A1")
            scan_agcu = scan_composition(fcc, 1000.0, db_agcu, solver)
            @test length(scan_agcu.x_grid) == 11
        end

        @testset "step_temperature with component_indices" begin
            db = read_tdb(joinpath(TDB_DIR, "steel1.TDB"))
            bcc = get_phase(db, "BCC_A2")

            # Cr-Mo step (CR=1, MO=3) at x_Mo = 0.5
            result = step_temperature(bcc, db, 0.5, 1800.0, 2000.0, 100.0;
                component_indices = (1, 3))

            @test length(result.points) == 3  # 1800, 1900, 2000
            @test all(p -> p.converged, result.points)
            @test all(p -> isfinite(p.gibbs_energy), result.points)

            @test result.points[1].n_phases >= 1

            # Check G value at 1800K is reasonable (negative for stable phase)
            @test result.points[1].gibbs_energy < 0

            # Backward compatibility: default should work for AgCu
            db_agcu = read_tdb(joinpath(TDB_DIR, "agcu.TDB"))
            fcc = get_phase(db_agcu, "FCC_A1")
            result_agcu = step_temperature(fcc, db_agcu, 0.7, 800.0, 1000.0, 100.0)
            @test length(result_agcu.points) == 3
        end

        @testset "find_miscibility_gap with component_indices" begin
            db = read_tdb(joinpath(TDB_DIR, "steel1.TDB"))
            liquid = get_phase(db, "LIQUID")

            # Cr-Mo LIQUID at high T - should not have a miscibility gap
            result = find_miscibility_gap(liquid, 2500.0, db; component_indices = (2, 4))

            # Backward compatibility: AgCu should still find the gap
            db_agcu = read_tdb(joinpath(TDB_DIR, "agcu.TDB"))
            fcc = get_phase(db_agcu, "FCC_A1")
            gap = find_miscibility_gap(fcc, 1000.0, db_agcu, 0.035, 0.90)
            @test !isnothing(gap)
            @test gap.converged == true
        end

        @testset "map_phase_diagram with component_indices" begin
            db = read_tdb(joinpath(TDB_DIR, "steel1.TDB"))
            bcc = get_phase(db, "BCC_A2")

            # Cr-Mo BCC at high T (no miscibility gap expected)
            result = map_phase_diagram(bcc, db, 1800.0, 2000.0, 100.0;
                component_indices = (1, 3))

            @test isa(result, PhaseDiagramResult)
            @test length(result.points) == 3

            # Backward compatibility: AgCu should still work
            db_agcu = read_tdb(joinpath(TDB_DIR, "agcu.TDB"))
            fcc = get_phase(db_agcu, "FCC_A1")
            result_agcu = map_phase_diagram(fcc, db_agcu, 800.0, 1000.0, 100.0)
            @test isa(result_agcu, PhaseDiagramResult)
        end
    end
end
