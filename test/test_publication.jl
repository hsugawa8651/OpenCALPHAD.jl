# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - PythonPlot extension tests

using Test
using OpenCALPHAD
using PythonPlot

# Force non-interactive backend on CI
PythonPlot.matplotlib.use("Agg")

@testset "savefig_publication (PythonPlot)" begin
    # ── Fixtures ──
    tdb_path = joinpath(dirname(@__DIR__), "reftest", "tdb", "agcu.TDB")
    db = read_tdb(tdb_path)
    fcc = get_phase(db, "FCC_A1")

    # GridScanResult: 11-point Gibbs curve at T=1000K
    grid_solver = GridSearchSolver(n_points = 11)
    grid_result = scan_composition(fcc, 1000.0, db, grid_solver)

    # StepResult: 3-point T scan with manual StepPoints
    sp1 = StepPoint(800.0,  -39427.0, 2, [0.275, 0.725],
                    [[0.009, 0.991], [0.962, 0.038]], true)
    sp2 = StepPoint(900.0,  -46335.0, 2, [0.256, 0.744],
                    [[0.018, 0.982], [0.935, 0.065]], true)
    sp3 = StepPoint(1000.0, -53648.0, 2, [0.228, 0.772],
                    [[0.034, 0.966], [0.897, 0.103]], true)
    step_result = StepResult(:T, [800.0, 900.0, 1000.0], [sp1, sp2, sp3], 0.7)

    # PhaseDiagramResult: 3-point binary
    pd1 = PhaseDiagramPoint(800.0,  (0.009, 0.962), true)
    pd2 = PhaseDiagramPoint(900.0,  (0.018, 0.935), true)
    pd3 = PhaseDiagramPoint(1000.0, (0.034, 0.897), true)
    pd_result = PhaseDiagramResult(fcc, [pd1, pd2, pd3])

    @testset "GridScanResult: Gibbs curve" begin
        mktempdir() do tmp
            path = joinpath(tmp, "gibbs.pdf")
            ret = savefig_publication(grid_result, path)
            @test isfile(path)
            @test ret == path
            @test filesize(path) > 0
        end
    end

    @testset "StepResult kind=:gibbs (default)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "step_g.pdf")
            savefig_publication(step_result, path)
            @test isfile(path)
        end
    end

    @testset "StepResult kind=:phase_fraction" begin
        mktempdir() do tmp
            path = joinpath(tmp, "step_pf.pdf")
            savefig_publication(step_result, path; kind = :phase_fraction)
            @test isfile(path)
        end
    end

    @testset "PhaseDiagramResult" begin
        mktempdir() do tmp
            path = joinpath(tmp, "pd.pdf")
            savefig_publication(pd_result, path)
            @test isfile(path)
        end
    end

    @testset "axis_width_cm / axis_height_cm custom kwargs" begin
        mktempdir() do tmp
            path = joinpath(tmp, "custom.pdf")
            savefig_publication(grid_result, path;
                                axis_width_cm = 10.0, axis_height_cm = 5.0)
            @test isfile(path)
        end
    end

    @testset "title kwarg" begin
        mktempdir() do tmp
            path = joinpath(tmp, "titled.pdf")
            savefig_publication(grid_result, path; title = "Ag-Cu @ 1000K")
            @test isfile(path)
        end
    end

    @testset "PNG output (extension dispatch)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "gibbs.png")
            savefig_publication(grid_result, path)
            @test isfile(path)
        end
    end

    @testset "kind unknown value throws ArgumentError" begin
        @test_throws ArgumentError savefig_publication(
            step_result, joinpath(tempdir(), "x.pdf"); kind = :bogus)
    end
end
