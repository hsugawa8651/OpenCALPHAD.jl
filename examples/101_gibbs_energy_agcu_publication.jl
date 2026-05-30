# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Ag-Cu Gibbs energy publication-quality figures
#
# Generates publication-quality PDF / PNG of the Ag-Cu FCC Gibbs energy curve
# at T = 1000 K, using the PythonPlot extension (matplotlib backend).
#
# Reuses the same scan_composition pattern as 101_gibbs_energy_agcu.jl, but
# emits PDF / PNG via savefig_publication for paper / poster figures.
#
# Usage:
#   julia --project=. examples/101_gibbs_energy_agcu_publication.jl
#
# Output (gitignored):
#   examples/101_gibbs_energy_agcu_publication.pdf
#   examples/101_gibbs_energy_agcu_publication.png

using OpenCALPHAD
using PythonPlot

# Force a non-interactive backend so this script works under CI / SSH
PythonPlot.matplotlib.use("Agg")

# Load Ag-Cu thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

# Publication-grade composition grid (51 points)
T = 1000.0  # K
solver = GridSearchSolver(n_points = 51)
result = scan_composition(fcc, T, db, solver)

stem = joinpath(@__DIR__, "101_gibbs_energy_agcu_publication")

savefig_publication(result, "$stem.pdf";
                    axis_width_cm = 8.0, axis_height_cm = 6.0,
                    title = "Ag-Cu FCC Gibbs energy at T = $(Int(T)) K",
                    show_minimum = true)
savefig_publication(result, "$stem.png";
                    axis_width_cm = 8.0, axis_height_cm = 6.0,
                    title = "Ag-Cu FCC Gibbs energy at T = $(Int(T)) K",
                    show_minimum = true)

println("Wrote $stem.pdf and $stem.png")
