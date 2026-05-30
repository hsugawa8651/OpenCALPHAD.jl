# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Documentation build script

using Documenter
using OpenCALPHAD
using PythonPlot                   # loads OpenCALPHADPythonPlotExt so its docstrings render
PythonPlot.matplotlib.use("Agg")  # headless backend for CI / @example figures

# The PythonPlot extension carries the docstrings for plot_on_axis! /
# figure_publication / savefig_publication; include it so @docs can render them.
const OCPythonPlotExt = Base.get_extension(OpenCALPHAD, :OpenCALPHADPythonPlotExt)

makedocs(
    sitename = "OpenCALPHAD.jl",
    authors = "Hiroharu Sugawara <hsugawa@gmail.com>",
    repo = Documenter.Remotes.GitHub("hsugawa8651", "OpenCALPHAD.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    modules = [OpenCALPHAD, OCPythonPlotExt],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Workflow" => "workflow.md",
        "Examples" => "examples.md",
        "Plotting" => "plotting.md",
        "Julia DSL" => "julia_dsl.md",
        "API Reference" => "api.md",
    ],
)

# Uncomment when deploying to GitHub Pages
deploydocs(
    repo = "github.com/hsugawa8651/OpenCALPHAD.jl.git",
    devbranch = "main",
)
