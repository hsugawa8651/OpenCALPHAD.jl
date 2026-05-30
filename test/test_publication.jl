# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - PythonPlot extension tests

using Test
using OpenCALPHAD
using PythonPlot

# Force non-interactive backend on CI
PythonPlot.matplotlib.use("Agg")

# Read a matplotlib title back as a Julia string (robust to Py wrapping)
_title(ax) = string(ax.get_title())

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

    @testset "axis_width_mm / axis_height_mm custom kwargs" begin
        mktempdir() do tmp
            path = joinpath(tmp, "custom.pdf")
            savefig_publication(grid_result, path;
                                axis_width_mm = 100.0, axis_height_mm = 50.0)
            @test isfile(path)
        end
    end

    @testset "axis_*_cm deprecated alias still works" begin
        mktempdir() do tmp
            path = joinpath(tmp, "cm.pdf")
            savefig_publication(grid_result, path;
                                axis_width_cm = 10.0, axis_height_cm = 5.0)
            @test isfile(path)
        end
    end

    @testset "mm takes precedence when both mm and cm given" begin
        # mm wins; the cm alias is ignored (with a deprecation warning), no error
        fig, ax = OpenCALPHAD.figure_publication(
            grid_result; axis_width_mm = 80.0, axis_width_cm = 8.0)
        @test fig isa PythonPlot.Figure
        PythonPlot.close(fig)
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
        # After the 3-layer refactor this ArgumentError originates in
        # plot_on_axis!(::StepResult); behavior unchanged for savefig callers.
        @test_throws ArgumentError savefig_publication(
            step_result, joinpath(tempdir(), "x.pdf"); kind = :bogus)
    end

    # ─────────────────── L1: plot_on_axis! (subplot composition) ───────────────────
    @testset "L1 plot_on_axis! returns ax (4 result paths)" begin
        for (r, kw) in ((grid_result, (;)),
                        (step_result, (;)),
                        (step_result, (kind = :phase_fraction,)),
                        (pd_result, (;)))
            fig = PythonPlot.figure()
            ax = fig.add_subplot()
            @test OpenCALPHAD.plot_on_axis!(ax, r; kw...) === ax
            PythonPlot.close(fig)
        end
    end

    @testset "L1 2x2 subplot composition smoke" begin
        mktempdir() do tmp
            fig = PythonPlot.figure(figsize = (10, 8))
            OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 1), grid_result)
            OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 2), step_result; kind = :gibbs)
            OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 3), step_result; kind = :phase_fraction)
            OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 4), pd_result)
            path = joinpath(tmp, "composite.pdf")
            fig.savefig(path)
            PythonPlot.close(fig)
            @test isfile(path)
        end
    end

    @testset "L1 plot_on_axis!(::StepResult; kind=:bogus) throws" begin
        fig = PythonPlot.figure()
        ax = fig.add_subplot()
        @test_throws ArgumentError OpenCALPHAD.plot_on_axis!(ax, step_result; kind = :bogus)
        PythonPlot.close(fig)
    end

    @testset "L1 title resolution (2-stage: empty -> no title, kwarg -> override)" begin
        fig = PythonPlot.figure(); ax = fig.add_subplot()
        OpenCALPHAD.plot_on_axis!(ax, grid_result)
        @test isempty(_title(ax))                       # no PF-style default title
        PythonPlot.close(fig)
        fig = PythonPlot.figure(); ax = fig.add_subplot()
        OpenCALPHAD.plot_on_axis!(ax, grid_result; title = "Ag-Cu")
        @test occursin("Ag-Cu", _title(ax))
        PythonPlot.close(fig)
    end

    # ─────────────────── L2: figure_publication ───────────────────
    @testset "L2 figure_publication returns (fig, ax) for all 3 result types" begin
        for r in (grid_result, step_result, pd_result)
            fig, ax = OpenCALPHAD.figure_publication(r)
            @test fig isa PythonPlot.Figure
            PythonPlot.close(fig)
        end
    end

    @testset "L2 StepResult kind flows through to L1" begin
        fig, ax = OpenCALPHAD.figure_publication(step_result; kind = :phase_fraction)
        @test fig isa PythonPlot.Figure
        PythonPlot.close(fig)
        @test_throws ArgumentError OpenCALPHAD.figure_publication(step_result; kind = :bogus)
    end

    @testset "L2 ylims routing" begin
        fig, ax = OpenCALPHAD.figure_publication(grid_result; ylims = (-60000.0, 0.0))
        @test fig isa PythonPlot.Figure
        @test occursin("-60000", string(ax.get_ylim()))
        PythonPlot.close(fig)
    end

    # ─────────────────── fallback (bad type → src ArgumentError) ───────────────────
    @testset "fallback: unsupported type throws ArgumentError" begin
        @test_throws ArgumentError OpenCALPHAD.plot_on_axis!(nothing, "not a result")
        @test_throws ArgumentError OpenCALPHAD.figure_publication("not a result")
    end
end
