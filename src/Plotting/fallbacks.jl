# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Fallback errors for extension-dependent functions

function savefig_publication(args...; kwargs...)
    throw(ArgumentError(
        "savefig_publication requires PythonPlot.jl. Run `using PythonPlot` first."
    ))
end
