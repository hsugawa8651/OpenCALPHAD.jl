# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Plot recipes for result types
#
# Recipe definitions moved to ext/OpenCALPHADRecipesBaseExt.jl
# Convenience functions implemented in ext/OpenCALPHADPlotsExt.jl

"""
    plot_phase_diagram(result::PhaseDiagramResult)

Plot a binary phase diagram. Requires `using Plots`.
"""
function plot_phase_diagram end

"""
    plot_gibbs_curve(result::GridScanResult)

Plot a Gibbs energy curve. Requires `using Plots`.
"""
function plot_gibbs_curve end

# PythonPlot extension API (3-layer). Implemented — with their authoritative,
# per-method docstrings — in ext/OpenCALPHADPythonPlotExt.jl. The extension
# module is added to Documenter's `modules` so `@docs` renders those docstrings.
function savefig_publication end
function plot_on_axis! end
function figure_publication end
