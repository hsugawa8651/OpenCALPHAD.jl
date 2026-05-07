# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - PythonPlot extension
#
# Convention: do NOT use `PythonPlot.subplots()` — it returns a raw
# `PythonCall.Py` whose `PythonPlot.close(fig)` triggers a `MethodError`.
# Always create figures via `PythonPlot.figure(figsize=...)` and add axes
# with `fig.add_axes([l, b, w, h])`.

module OpenCALPHADPythonPlotExt

using OpenCALPHAD
import OpenCALPHAD: GridScanResult, StepResult, PhaseDiagramResult,
                    temperatures, gibbs_energies, phase_fractions,
                    savefig_publication
import PythonPlot

const CM_PER_INCH = 2.54

# axis_variable -> default xlabel mapping for StepResult
const AXIS_LABEL = Dict(
    :T => "Temperature T [K]",
    :P => "Pressure P [Pa]",
    :x => "Composition x",
)

"""
    _layout_axes(axis_width_cm, axis_height_cm, n;
                 margin_left_cm = 3.0, margin_right_cm = 0.3,
                 margin_bottom_cm = 1.5, margin_top_cm = 0.8,
                 hgap_cm = 1.8, vgap_cm = 1.5,
                 nrows = 1, ncols = 1)

Compute figure size (inches) and axes positions from axis dimensions and layout.

`margin_left_cm = 3.0` (widened from BT/PX default 2.0) accommodates OC's 4-5 digit
G [J/mol] tick labels plus the y-axis label.
"""
function _layout_axes(
    axis_width_cm,
    axis_height_cm,
    n;
    margin_left_cm   = 3.0,
    margin_right_cm  = 0.3,
    margin_bottom_cm = 1.5,
    margin_top_cm    = 0.8,
    hgap_cm          = 1.8,
    vgap_cm          = 1.5,
    nrows            = 1,
    ncols            = 1,
)
    widths  = axis_width_cm  isa AbstractVector ? axis_width_cm  : fill(axis_width_cm,  ncols)
    heights = axis_height_cm isa AbstractVector ? axis_height_cm : fill(axis_height_cm, nrows)

    fig_w_cm = margin_left_cm   + sum(widths)  + hgap_cm * (ncols - 1) + margin_right_cm
    fig_h_cm = margin_bottom_cm + sum(heights) + vgap_cm * (nrows - 1) + margin_top_cm

    fig_w = fig_w_cm / CM_PER_INCH
    fig_h = fig_h_cm / CM_PER_INCH

    positions = Vector{NTuple{4,Float64}}()
    for row in 1:nrows
        for col in 1:ncols
            left   = (margin_left_cm   + sum(widths[1:(col - 1)])    + hgap_cm * (col - 1))     / fig_w_cm
            bottom = (margin_bottom_cm + sum(heights[(row + 1):end]) + vgap_cm * (nrows - row)) / fig_h_cm
            w = widths[col]  / fig_w_cm
            h = heights[row] / fig_h_cm
            push!(positions, (left, bottom, w, h))
        end
    end
    return fig_w, fig_h, positions
end

"""
    _plot_gibbs_curve_on_ax!(ax, r::GridScanResult; kwargs...)

Plot a `GridScanResult` Gibbs energy curve on the given matplotlib `ax`.
Empty-string `xlabel` / `ylabel` / `title` use defaults; non-empty values override.
`show_minimum = true` draws a vertical line at the minimum-G composition.
"""
function _plot_gibbs_curve_on_ax!(
    ax,
    r::GridScanResult;
    xlabel::AbstractString = "",
    ylabel::AbstractString = "",
    title::AbstractString  = "",
    color = "black",
    linewidth = 2.0,
    linestyle = "-",
    show_minimum::Bool = false,
)
    ax.plot(r.x_grid, r.G_values; color = color, linewidth = linewidth, linestyle = linestyle)
    if show_minimum
        min_idx = argmin(r.G_values)
        ax.axvline(r.x_grid[min_idx]; color = "gray", linewidth = 0.5, linestyle = "--")
    end
    eff_xlabel = isempty(xlabel) ? "Composition x" : xlabel
    eff_ylabel = isempty(ylabel) ? "G [J/mol]"     : ylabel
    ax.set_xlabel(eff_xlabel)
    ax.set_ylabel(eff_ylabel)
    isempty(title) || ax.set_title(title)
    return ax
end

"""
    _plot_step_gibbs_on_ax!(ax, r::StepResult; kwargs...)

Plot a `StepResult` Gibbs energy curve along its axis variable.
The default xlabel uses `AXIS_LABEL[r.axis_variable]`.
"""
function _plot_step_gibbs_on_ax!(
    ax,
    r::StepResult;
    xlabel::AbstractString = "",
    ylabel::AbstractString = "",
    title::AbstractString  = "",
    color = "black",
    linewidth = 2.0,
    linestyle = "-",
)
    ax.plot(temperatures(r), gibbs_energies(r);
            color = color, linewidth = linewidth, linestyle = linestyle)
    eff_xlabel = isempty(xlabel) ?
        get(AXIS_LABEL, r.axis_variable, String(r.axis_variable)) : xlabel
    eff_ylabel = isempty(ylabel) ? "G [J/mol]" : ylabel
    ax.set_xlabel(eff_xlabel)
    ax.set_ylabel(eff_ylabel)
    isempty(title) || ax.set_title(title)
    return ax
end

"""
    _plot_step_phase_fraction_on_ax!(ax, r::StepResult; kwargs...)

Plot two phase fractions of a `StepResult` along its axis variable.
Forces `ylim = (0, 1)` and adds a legend at `legend_loc`.
"""
function _plot_step_phase_fraction_on_ax!(
    ax,
    r::StepResult;
    xlabel::AbstractString = "",
    ylabel::AbstractString = "",
    title::AbstractString  = "",
    color1 = "C0",
    color2 = "C1",
    linewidth = 2.0,
    legend_loc::AbstractString = "center right",
)
    xs = temperatures(r)
    f1 = phase_fractions(r, 1)
    f2 = phase_fractions(r, 2)
    ax.plot(xs, f1; color = color1, linewidth = linewidth, label = "phase 1")
    ax.plot(xs, f2; color = color2, linewidth = linewidth, label = "phase 2")
    ax.set_ylim(0, 1)
    ax.legend(loc = legend_loc)
    eff_xlabel = isempty(xlabel) ?
        get(AXIS_LABEL, r.axis_variable, String(r.axis_variable)) : xlabel
    eff_ylabel = isempty(ylabel) ? "Phase fraction" : ylabel
    ax.set_xlabel(eff_xlabel)
    ax.set_ylabel(eff_ylabel)
    isempty(title) || ax.set_title(title)
    return ax
end

"""
    _plot_phase_diagram_on_ax!(ax, r::PhaseDiagramResult; kwargs...)

Plot a binary phase diagram from a `PhaseDiagramResult`. The two-phase region
between the boundaries is shaded with `fill_color` at `fill_alpha`. Skips
non-converged points; if no points converge, an empty axes is drawn with axis
labels only.
"""
function _plot_phase_diagram_on_ax!(
    ax,
    r::PhaseDiagramResult;
    xlabel::AbstractString = "",
    ylabel::AbstractString = "",
    title::AbstractString  = "",
    boundary_color = "black",
    linewidth = 2.0,
    fill_color::AbstractString = "lightgray",
    fill_alpha::Real = 0.3,
)
    eff_xlabel = isempty(xlabel) ? "Composition x" : xlabel
    eff_ylabel = isempty(ylabel) ? "T [K]"        : ylabel

    converged_pts = filter(p -> p.converged, r.points)
    if isempty(converged_pts)
        ax.set_xlabel(eff_xlabel)
        ax.set_ylabel(eff_ylabel)
        isempty(title) || ax.set_title(title)
        return ax
    end

    T_vals  = [p.temperature     for p in converged_pts]
    x1_vals = [p.compositions[1] for p in converged_pts]
    x2_vals = [p.compositions[2] for p in converged_pts]

    ax.fill_betweenx(T_vals, x1_vals, x2_vals; alpha = fill_alpha, color = fill_color)
    ax.plot(x1_vals, T_vals; color = boundary_color, linewidth = linewidth)
    ax.plot(x2_vals, T_vals; color = boundary_color, linewidth = linewidth)

    ax.set_xlabel(eff_xlabel)
    ax.set_ylabel(eff_ylabel)
    isempty(title) || ax.set_title(title)
    return ax
end

# ── savefig_publication dispatch ──

"""
    OpenCALPHAD.savefig_publication(r::GridScanResult, path; kwargs...) -> path

Render `r` as a publication-quality figure to `path` (PDF/PNG; format chosen by
extension). Returns `path` for chaining.
"""
function OpenCALPHAD.savefig_publication(
    r::GridScanResult,
    path::AbstractString;
    axis_width_cm = 8.0,
    axis_height_cm = 6.0,
    title::AbstractString = "",
    ylims = nothing,
    kwargs...,
)
    fig_w, fig_h, positions = _layout_axes(axis_width_cm, axis_height_cm, 1)
    fig = PythonPlot.figure(figsize = (fig_w, fig_h))
    ax  = fig.add_axes(collect(positions[1]))
    _plot_gibbs_curve_on_ax!(ax, r; title = title, kwargs...)
    isnothing(ylims) || ax.set_ylim(ylims...)
    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

"""
    OpenCALPHAD.savefig_publication(r::StepResult, path; kind = :gibbs, kwargs...) -> path

Render a `StepResult` to `path`. `kind = :gibbs` plots `G(axis)`, `kind = :phase_fraction`
plots phase fractions vs. axis. Other `kind` values raise `ArgumentError`.
"""
function OpenCALPHAD.savefig_publication(
    r::StepResult,
    path::AbstractString;
    kind::Symbol = :gibbs,
    axis_width_cm = 8.0,
    axis_height_cm = 6.0,
    title::AbstractString = "",
    ylims = nothing,
    kwargs...,
)
    fig_w, fig_h, positions = _layout_axes(axis_width_cm, axis_height_cm, 1)
    fig = PythonPlot.figure(figsize = (fig_w, fig_h))
    ax  = fig.add_axes(collect(positions[1]))
    if kind === :gibbs
        _plot_step_gibbs_on_ax!(ax, r; title = title, kwargs...)
    elseif kind === :phase_fraction
        _plot_step_phase_fraction_on_ax!(ax, r; title = title, kwargs...)
    else
        PythonPlot.close(fig)
        throw(ArgumentError("unknown kind=$(kind); expected :gibbs or :phase_fraction"))
    end
    isnothing(ylims) || ax.set_ylim(ylims...)
    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

"""
    OpenCALPHAD.savefig_publication(r::PhaseDiagramResult, path; kwargs...) -> path

Render a `PhaseDiagramResult` (binary phase diagram with shaded two-phase region).
"""
function OpenCALPHAD.savefig_publication(
    r::PhaseDiagramResult,
    path::AbstractString;
    axis_width_cm = 8.0,
    axis_height_cm = 6.0,
    title::AbstractString = "",
    ylims = nothing,
    kwargs...,
)
    fig_w, fig_h, positions = _layout_axes(axis_width_cm, axis_height_cm, 1)
    fig = PythonPlot.figure(figsize = (fig_w, fig_h))
    ax  = fig.add_axes(collect(positions[1]))
    _plot_phase_diagram_on_ax!(ax, r; title = title, kwargs...)
    isnothing(ylims) || ax.set_ylim(ylims...)
    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

end # module
