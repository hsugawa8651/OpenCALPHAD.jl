# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: composite publication figure via plot_on_axis!
#
# Demonstrates the 3-layer PythonPlot API by composing four Ag-Cu result types
# into a single 2x2 publication panel using OpenCALPHAD.plot_on_axis! (L1):
#   (1) GridScanResult     - Gibbs energy curve G(x) at 1000 K
#   (2) StepResult         - G(T) along a STEP calculation (kind = :gibbs)
#   (3) StepResult         - phase fractions vs T   (kind = :phase_fraction)
#   (4) PhaseDiagramResult - binary phase diagram (FCC miscibility gap)
#
# The caller owns the figure (created with PythonPlot.figure / add_subplot) and
# is responsible for PythonPlot.close. Candidate figure for the JuliaCon 2026
# Mainz poster (AFVQLM).
#
# Usage:
#   julia --project=. examples/811_composite_publication.jl
#
# Output (gitignored):
#   examples/811_composite_publication.pdf
#   examples/811_composite_publication.png

using OpenCALPHAD
using PythonPlot

PythonPlot.matplotlib.use("Agg")

# Larger fonts for poster / publication readability
let rc = PythonPlot.matplotlib.rcParams
    rc["axes.titlesize"]  = 17
    rc["axes.labelsize"]  = 16
    rc["xtick.labelsize"] = 14
    rc["ytick.labelsize"] = 14
    rc["legend.fontsize"] = 14
end

tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "agcu.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

# Four result types (same Ag-Cu system / ranges as the baseline examples)
grid_result = scan_composition(fcc, 1000.0, db, GridSearchSolver(n_points = 51))
step_result = step_temperature(fcc, db, 0.5, 800.0, 1300.0, 25.0)
pd_result   = map_phase_diagram(fcc, db, 800.0, 1300.0, 25.0)

# Compose a 2x2 grid: caller creates the figure and adds the axes (L1 only
# mutates each axis). subplots() on the caller side is fine — we own `close`.
fig = PythonPlot.figure(figsize = (12, 9))
OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 1), grid_result;
                          title = "Gibbs energy at 1000 K", show_minimum = true)
OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 2), step_result;
                          kind = :gibbs, title = "STEP: G vs T")
OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 3), step_result;
                          kind = :phase_fraction, title = "Phase fractions vs T",
                          legend_loc = "lower center")
OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 4), pd_result;
                          title = "Binary phase diagram")
fig.tight_layout()

stem = joinpath(@__DIR__, "811_composite_publication")
fig.savefig("$stem.pdf")
fig.savefig("$stem.png")
PythonPlot.close(fig)

println("Wrote $stem.pdf and $stem.png")
