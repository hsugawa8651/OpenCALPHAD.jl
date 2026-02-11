# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Tests for multi-component support

@testset "Multi-component support" begin
    @testset "find_va_index" begin
        # VA is first element
        @test OpenCALPHAD.find_va_index(["VA"]) == 1

        # VA is second element (steel1.TDB case)
        @test OpenCALPHAD.find_va_index(["C", "VA"]) == 2

        # VA is third element
        @test OpenCALPHAD.find_va_index(["A", "B", "VA"]) == 3

        # No VA - fallback to 1
        @test OpenCALPHAD.find_va_index(["CR", "FE"]) == 1
    end

    @testset "make_y_matrix" begin
        # Test with steel1.TDB (Cr-Mo)
        db = read_tdb(joinpath(@__DIR__, "..", "reftest", "tdb", "steel1.TDB"))
        bcc = get_phase(db, "BCC_A2")

        # Cr-Mo (indices 1, 3) at x_Mo = 0.5
        y = OpenCALPHAD.make_y_matrix(bcc, 0.5, (1, 3))
        @test y[1, 1] ≈ 0.5  # Cr (1 - x)
        @test y[1, 3] ≈ 0.5  # Mo (x)
        @test y[2, 2] ≈ 1.0  # VA (index 2 in sublattice 2: ["C", "VA"])

        # Test with agcu.TDB (Ag-Cu) - default indices (1, 2)
        db_agcu = read_tdb(joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB"))
        fcc = get_phase(db_agcu, "FCC_A1")

        y2 = OpenCALPHAD.make_y_matrix(fcc, 0.3, (1, 2))
        @test y2[1, 1] ≈ 0.7  # Ag (1 - x)
        @test y2[1, 2] ≈ 0.3  # Cu (x)
        @test y2[2, 1] ≈ 1.0  # VA (index 1 in sublattice 2: ["VA"])

        # Test default component_indices
        y3 = OpenCALPHAD.make_y_matrix(fcc, 0.3)
        @test y3[1, 1] ≈ 0.7
        @test y3[1, 2] ≈ 0.3
    end
end
