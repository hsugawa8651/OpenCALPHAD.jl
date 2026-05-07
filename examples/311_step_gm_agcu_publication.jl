# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Ag-Cu STEP publication-quality figures
#
# Generates publication-quality PDF / PNG of an Ag-Cu STEP calculation, both
# kinds of plot (Gibbs energy and phase fractions vs T), via the PythonPlot
# extension (matplotlib backend).
#
# Note: this script uses `step_temperature` (returning a `StepResult`) to
# demonstrate `savefig_publication(::StepResult)` dispatch with `kind = :gibbs`
# and `kind = :phase_fraction`. The companion baseline 311_step_gm_agcu.jl
# uses `scan_composition` + `thermodynamic_properties` for a 4-panel Plots.jl
# figure (different shape — same Ag-Cu STEP topic).
#
# Usage:
#   julia --project=. examples/311_step_gm_agcu_publication.jl
#
# Output (gitignored):
#   examples/311_step_gm_agcu_gibbs_publication.pdf
#   examples/311_step_gm_agcu_gibbs_publication.png
#   examples/311_step_gm_agcu_phase_fraction_publication.pdf
#   examples/311_step_gm_agcu_phase_fraction_publication.png

using OpenCALPHAD
using PythonPlot

PythonPlot.matplotlib.use("Agg")

tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

# STEP over T = 800..1300 K, x_overall = 0.5 (Ag-Cu midpoint, well inside the gap)
step_result = step_temperature(fcc, db, 0.5, 800.0, 1300.0, 25.0)

stem = joinpath(@__DIR__, "311_step_gm_agcu")

for kind in (:gibbs, :phase_fraction)
    for ext in ("pdf", "png")
        savefig_publication(step_result, "$(stem)_$(kind)_publication.$(ext)";
                            kind = kind,
                            axis_width_cm = 8.0, axis_height_cm = 6.0,
                            title = "Ag-Cu STEP, x(Cu) = 0.5")
    end
end

println("Wrote 4 files under $(stem)_*_publication.{pdf,png}")
