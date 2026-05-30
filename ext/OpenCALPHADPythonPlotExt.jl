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
                    savefig_publication, plot_on_axis!, figure_publication
import PythonPlot

const MM_PER_INCH = 25.4

# Resolve an axis dimension in mm from the mm keyword (canonical) and the
# deprecated cm alias. The mm keyword always takes precedence; the cm alias
# emits a deprecation warning, and is ignored when the mm form is also given.
function _resolve_dim_mm(mm, cm, default, name)
    if cm !== nothing
        if mm !== nothing
            Base.depwarn("axis_$(name)_cm is deprecated and ignored because " *
                         "axis_$(name)_mm was also given; use axis_$(name)_mm only",
                         :figure_publication)
            return mm
        end
        Base.depwarn("axis_$(name)_cm is deprecated; use axis_$(name)_mm (mm = 10 × cm)",
                     :figure_publication)
        return 10.0 * cm
    end
    return mm === nothing ? default : mm
end

# axis_variable -> default xlabel mapping for StepResult
const AXIS_LABEL = Dict(
    :T => "Temperature T [K]",
    :P => "Pressure P [Pa]",
    :x => "Composition x",
)

"""
    _layout_axes(axis_width_mm, axis_height_mm, n;
                 margin_left_mm = 30.0, margin_right_mm = 3.0,
                 margin_bottom_mm = 15.0, margin_top_mm = 8.0,
                 hgap_mm = 18.0, vgap_mm = 15.0,
                 nrows = 1, ncols = 1)

Compute figure size (inches) and axes positions from axis dimensions (in mm)
and layout.

`margin_left_mm = 30.0` (widened from BT/PX default 20.0) accommodates OC's
4-5 digit G [J/mol] tick labels plus the y-axis label.
"""
function _layout_axes(
    axis_width_mm,
    axis_height_mm,
    n;
    margin_left_mm   = 30.0,
    margin_right_mm  = 3.0,
    margin_bottom_mm = 15.0,
    margin_top_mm    = 8.0,
    hgap_mm          = 18.0,
    vgap_mm          = 15.0,
    nrows            = 1,
    ncols            = 1,
)
    widths  = axis_width_mm  isa AbstractVector ? axis_width_mm  : fill(axis_width_mm,  ncols)
    heights = axis_height_mm isa AbstractVector ? axis_height_mm : fill(axis_height_mm, nrows)

    fig_w_mm = margin_left_mm   + sum(widths)  + hgap_mm * (ncols - 1) + margin_right_mm
    fig_h_mm = margin_bottom_mm + sum(heights) + vgap_mm * (nrows - 1) + margin_top_mm

    fig_w = fig_w_mm / MM_PER_INCH
    fig_h = fig_h_mm / MM_PER_INCH

    positions = Vector{NTuple{4,Float64}}()
    for row in 1:nrows
        for col in 1:ncols
            left   = (margin_left_mm   + sum(widths[1:(col - 1)])    + hgap_mm * (col - 1))     / fig_w_mm
            bottom = (margin_bottom_mm + sum(heights[(row + 1):end]) + vgap_mm * (nrows - row)) / fig_h_mm
            w = widths[col]  / fig_w_mm
            h = heights[row] / fig_h_mm
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

# ══════════════════════════════════════════════════════════════════════
# L1 — plot_on_axis!(ax, r; ...) : mutate a user-supplied axis, return ax
# (thin public wrappers over the existing _plot_*_on_ax! helpers)
# ══════════════════════════════════════════════════════════════════════

"""
    plot_on_axis!(ax, r::GridScanResult; kwargs...) -> ax

Draw the Gibbs energy curve `G(x)` of a [`GridScanResult`](@ref) onto the
matplotlib axis you supply and return `ax`, so panels can be composed.
Requires `using PythonPlot`; not exported (call as `OpenCALPHAD.plot_on_axis!`).

# Arguments

- `ax`: a matplotlib axis (PythonCall `Py`) to draw onto
- `r::GridScanResult`: the composition scan to plot

# Keywords

- `color="black"`: curve color
- `linewidth=2.0`: curve line width
- `linestyle="-"`: curve line style
- `show_minimum=false`: draw a vertical line at the composition of minimum G
- `xlabel::AbstractString=""`: x axis label; empty uses `"Composition x"`
- `ylabel::AbstractString=""`: y axis label; empty uses `"G [J/mol]"`
- `title::AbstractString=""`: title; empty leaves the axes untitled
"""
OpenCALPHAD.plot_on_axis!(ax, r::GridScanResult; kwargs...) =
    _plot_gibbs_curve_on_ax!(ax, r; kwargs...)

"""
    plot_on_axis!(ax, r::PhaseDiagramResult; kwargs...) -> ax

Draw a binary phase diagram of a [`PhaseDiagramResult`](@ref) onto `ax` and
return `ax`; the region between the two boundaries is shaded. Points that did
not converge are skipped, and if none converge an empty labelled axes is drawn.
Requires `using PythonPlot`; not exported (call as `OpenCALPHAD.plot_on_axis!`).

# Arguments

- `ax`: a matplotlib axis (PythonCall `Py`) to draw onto
- `r::PhaseDiagramResult`: the mapped phase boundaries to plot

# Keywords

- `boundary_color="black"`: boundary line color
- `linewidth=2.0`: boundary line width
- `fill_color="lightgray"`: shading color of the region between the boundaries
- `fill_alpha=0.3`: shading opacity
- `xlabel::AbstractString=""`: x axis label; empty uses `"Composition x"`
- `ylabel::AbstractString=""`: y axis label; empty uses `"T [K]"`
- `title::AbstractString=""`: title; empty leaves the axes untitled
"""
OpenCALPHAD.plot_on_axis!(ax, r::PhaseDiagramResult; kwargs...) =
    _plot_phase_diagram_on_ax!(ax, r; kwargs...)

"""
    plot_on_axis!(ax, r::StepResult; kind=:gibbs, kwargs...) -> ax

Draw a [`StepResult`](@ref) onto `ax` and return `ax`. Requires
`using PythonPlot`; not exported (call as `OpenCALPHAD.plot_on_axis!`).

# Arguments

- `ax`: a matplotlib axis (PythonCall `Py`) to draw onto
- `r::StepResult`: the STEP calculation to plot

# Keywords

- `kind::Symbol=:gibbs`: `:gibbs` plots `G` against the axis variable;
    `:phase_fraction` plots the two phase fractions (forcing `ylim=(0,1)` with a
    legend). Any other value raises `ArgumentError`.
- `color`, `linewidth=2.0`, `linestyle="-"`: curve style when `kind=:gibbs`
- `color1="C0"`, `color2="C1"`: per phase colors when `kind=:phase_fraction`
- `legend_loc="center right"`: legend position when `kind=:phase_fraction`
- `xlabel::AbstractString=""`: x axis label; empty uses the axis variable label
- `ylabel::AbstractString=""`: y axis label; empty uses `"G [J/mol]"` or
    `"Phase fraction"`
- `title::AbstractString=""`: title; empty leaves the axes untitled
"""
function OpenCALPHAD.plot_on_axis!(ax, r::StepResult; kind::Symbol = :gibbs, kwargs...)
    if kind === :gibbs
        _plot_step_gibbs_on_ax!(ax, r; kwargs...)
    elseif kind === :phase_fraction
        _plot_step_phase_fraction_on_ax!(ax, r; kwargs...)
    else
        throw(ArgumentError("unknown kind=$(kind); expected :gibbs or :phase_fraction"))
    end
    return ax
end

# ══════════════════════════════════════════════════════════════════════
# L2 — figure_publication(r; ...) : create (fig, ax), single axis
# ══════════════════════════════════════════════════════════════════════

"""
    figure_publication(r; axis_width_mm=80.0, axis_height_mm=60.0, ylims=nothing, kwargs...) -> (fig, ax)

Create a publication matplotlib figure and a single axis for `r`, draw it via
[`plot_on_axis!`](@ref), and return `(fig, ax)` so you can tweak it before
saving. The caller owns the figure and must call `PythonPlot.close(fig)`.
Requires `using PythonPlot`; not exported (call as `OpenCALPHAD.figure_publication`).

# Arguments

- `r`: a `GridScanResult`, `StepResult`, or `PhaseDiagramResult`

# Keywords

- `axis_width_mm=80.0`: plotting area width in mm
- `axis_height_mm=60.0`: plotting area height in mm
- `axis_width_cm`, `axis_height_cm`: deprecated aliases (1 cm = 10 mm). The mm
    form takes precedence; the cm form warns and is ignored if both are given
- `ylims=nothing`: pass a tuple to override the y axis limits
- `kwargs...`: forwarded to [`plot_on_axis!`](@ref)
"""
function OpenCALPHAD.figure_publication(
    r::Union{GridScanResult,StepResult,PhaseDiagramResult};
    axis_width_mm = nothing,
    axis_height_mm = nothing,
    axis_width_cm = nothing,
    axis_height_cm = nothing,
    ylims = nothing,
    kwargs...,
)
    w_mm = _resolve_dim_mm(axis_width_mm, axis_width_cm, 80.0, "width")
    h_mm = _resolve_dim_mm(axis_height_mm, axis_height_cm, 60.0, "height")
    fig_w, fig_h, positions = _layout_axes(w_mm, h_mm, 1)
    fig = PythonPlot.figure(figsize = (fig_w, fig_h))
    try
        ax = fig.add_axes(collect(positions[1]))
        plot_on_axis!(ax, r; kwargs...)
        isnothing(ylims) || ax.set_ylim(ylims...)
        return fig, ax
    catch
        PythonPlot.close(fig)
        rethrow()
    end
end

# ══════════════════════════════════════════════════════════════════════
# L3 — savefig_publication(r, path; ...) : L2 + savefig + close
# ══════════════════════════════════════════════════════════════════════

"""
    savefig_publication(r, path; kwargs...) -> path

Render `r` to `path` and return `path`; the output format follows the file
extension (`.pdf` or `.png`). This is the convenience entry point (the only
exported layer): it delegates to [`figure_publication`](@ref) and closes the
figure for you. Requires `using PythonPlot`.

# Arguments

- `r`: a `GridScanResult`, `StepResult`, or `PhaseDiagramResult`
- `path::AbstractString`: output file; the extension selects PDF or PNG

# Keywords

- `kwargs...`: forwarded to [`figure_publication`](@ref) (and onward to
    [`plot_on_axis!`](@ref))
"""
function OpenCALPHAD.savefig_publication(
    r::Union{GridScanResult,StepResult,PhaseDiagramResult},
    path::AbstractString;
    kwargs...,
)
    fig, _ = figure_publication(r; kwargs...)
    try
        fig.savefig(path)
    finally
        PythonPlot.close(fig)
    end
    return path
end

end # module
