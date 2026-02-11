# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Documentation build script

using Documenter
using OpenCALPHAD

makedocs(
    sitename = "OpenCALPHAD.jl",
    authors = "Hiroharu Sugawara <hsugawa@gmail.com>",
    repo = Documenter.Remotes.GitHub("hsugawa8651", "OpenCALPHAD.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    modules = [OpenCALPHAD],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Workflow" => "workflow.md",
        "Examples" => "examples.md",
        "Julia DSL" => "julia_dsl.md",
        "API Reference" => "api.md",
    ],
)

# Uncomment when deploying to GitHub Pages
deploydocs(
    repo = "github.com/hsugawa8651/OpenCALPHAD.jl.git",
    devbranch = "main",
)
