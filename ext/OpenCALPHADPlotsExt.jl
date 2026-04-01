# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Plots.jl extension

module OpenCALPHADPlotsExt

using Plots
import OpenCALPHAD
using OpenCALPHAD: PhaseDiagramResult, GridScanResult

function OpenCALPHAD.plot_phase_diagram(result::PhaseDiagramResult)
    Plots.plot(result)
end

function OpenCALPHAD.plot_gibbs_curve(result::GridScanResult)
    Plots.plot(result)
end

end # module
