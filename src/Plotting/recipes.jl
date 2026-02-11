# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2012-2026 Bo Sundman and OpenCALPHAD contributors
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of OpenCALPHAD.jl - Plot recipes for result types

using RecipesBase

"""
    @recipe for GridScanResult

Plot Gibbs energy curve G(x) from grid scan result.
"""
@recipe function f(result::GridScanResult)
    xlabel --> "Composition x"
    ylabel --> "Gibbs Energy G [J/mol]"
    title --> "Gibbs Energy Curve"
    legend --> false
    linewidth --> 2

    result.x_grid, result.G_values
end

"""
    @recipe for StepResult

Plot Gibbs energy vs temperature from STEP result.
"""
@recipe function f(result::StepResult)
    xlabel --> "Temperature T [K]"
    ylabel --> "Gibbs Energy G [J/mol]"
    title --> "STEP Calculation"
    legend --> false
    linewidth --> 2

    temperatures(result), gibbs_energies(result)
end

"""
    @recipe for StepResult with :phase_fraction plottype

Plot phase fractions vs temperature.
Usage: plot(result, plottype=:phase_fraction)
"""
@recipe function f(result::StepResult, ::Val{:phase_fraction})
    xlabel --> "Temperature T [K]"
    ylabel --> "Phase Fraction"
    title --> "Phase Fractions vs Temperature"
    legend --> :topright
    linewidth --> 2
    ylims --> (0, 1)

    T = temperatures(result)
    f1 = phase_fractions(result, 1)
    f2 = phase_fractions(result, 2)

    @series begin
        label --> "Phase 1"
        T, f1
    end

    @series begin
        label --> "Phase 2"
        T, f2
    end
end

"""
    @recipe for PhaseDiagramResult

Plot binary phase diagram (T-x diagram) showing miscibility gap boundaries.
"""
@recipe function f(result::PhaseDiagramResult)
    xlabel --> "Composition x"
    ylabel --> "Temperature T [K]"
    title --> "Binary Phase Diagram"
    legend --> false
    linewidth --> 2

    # Filter converged points only
    converged = filter(p -> p.converged, result.points)

    if isempty(converged)
        # Return empty plot if no converged points
        return Float64[], Float64[]
    end

    # Extract data for left and right boundaries
    T_vals = [p.temperature for p in converged]
    x1_vals = [p.compositions[1] for p in converged]
    x2_vals = [p.compositions[2] for p in converged]

    # Create closed boundary curve
    # Go up left boundary, then down right boundary
    x_boundary = vcat(x1_vals, reverse(x2_vals))
    T_boundary = vcat(T_vals, reverse(T_vals))

    # Fill the two-phase region
    fillrange --> 0
    fillalpha --> 0.3
    seriestype --> :path

    x_boundary, T_boundary
end

"""
    plot_phase_diagram(result::PhaseDiagramResult; kwargs...)

Convenience function to plot a phase diagram with appropriate defaults.
Returns the plot recipe data for use with Plots.jl.
"""
function plot_phase_diagram(result::PhaseDiagramResult)
    return result
end

"""
    plot_gibbs_curve(result::GridScanResult; kwargs...)

Convenience function to plot Gibbs energy curve.
Returns the plot recipe data for use with Plots.jl.
"""
function plot_gibbs_curve(result::GridScanResult)
    return result
end
