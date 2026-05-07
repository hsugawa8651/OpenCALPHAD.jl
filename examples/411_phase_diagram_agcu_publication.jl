# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Ag-Cu phase diagram publication-quality figures
#
# Generates publication-quality PDF / PNG of the Ag-Cu binary phase diagram
# (FCC miscibility gap) via the PythonPlot extension (matplotlib backend).
#
# Reuses the same map_phase_diagram pattern as 411_phase_diagram_agcu.jl, but
# emits PDF / PNG via savefig_publication for paper / poster figures.
#
# Usage:
#   julia --project=. examples/411_phase_diagram_agcu_publication.jl
#
# Output (gitignored):
#   examples/411_phase_diagram_agcu_publication.pdf
#   examples/411_phase_diagram_agcu_publication.png

using OpenCALPHAD
using PythonPlot

PythonPlot.matplotlib.use("Agg")

tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

# Phase diagram from 800..1300 K (matches baseline 411 range)
result = map_phase_diagram(fcc, db, 800.0, 1300.0, 25.0)

stem = joinpath(@__DIR__, "411_phase_diagram_agcu_publication")

savefig_publication(result, "$stem.pdf";
                    axis_width_cm = 10.0, axis_height_cm = 7.0,
                    title = "Ag-Cu binary phase diagram (FCC miscibility gap)")
savefig_publication(result, "$stem.png";
                    axis_width_cm = 10.0, axis_height_cm = 7.0,
                    title = "Ag-Cu binary phase diagram (FCC miscibility gap)")

println("Wrote $stem.pdf and $stem.png")
