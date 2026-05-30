# Plotting

OpenCALPHAD.jl provides three plotting paths via package extensions:

| Extension | Trigger | Functionality |
|-----------|---------|---------------|
| RecipesBaseExt | `using RecipesBase` (or any backend) | `plot(result)` recipes for all result types |
| PlotsExt | `using Plots` | Convenience functions `plot_gibbs_curve`, `plot_phase_diagram` |
| PythonPlotExt | `using PythonPlot` | 3-layer API `plot_on_axis!` / `figure_publication` / `savefig_publication` for publication-quality static PDF / PNG |

---

## RecipesBase Recipes

These recipes are activated automatically when any RecipesBase-compatible backend is loaded.

```julia
using OpenCALPHAD
using Plots   # or CairoMakie, etc.
```

### Gibbs Energy Curve — `GridScanResult`

Plot Gibbs energy G(x) from a composition scan.

```julia
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")
solver = GridSearchSolver(n_points=101)
scan = scan_composition(fcc, 1000.0, db, solver)

plot(scan)
```

| Attribute | Default |
|-----------|---------|
| xlabel | `"Composition x"` |
| ylabel | `"Gibbs Energy G [J/mol]"` |
| title | `"Gibbs Energy Curve"` |
| linewidth | 2 |

Override any attribute with standard keyword arguments:

```julia
plot(scan, title="Ag-Cu FCC at 1000K", ylabel="G [kJ/mol]")
```

See also: [`examples/122_scan_crmo.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/122_scan_crmo.jl)

### Gibbs Energy vs Temperature — `StepResult`

Plot Gibbs energy G(T) from a STEP calculation.

```julia
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")
result = step_temperature(fcc, db, 0.5, 300.0, 1400.0, 10.0)

plot(result)
```

| Attribute | Default |
|-----------|---------|
| xlabel | `"Temperature T [K]"` |
| ylabel | `"Gibbs Energy G [J/mol]"` |
| title | `"STEP Calculation"` |
| linewidth | 2 |

See also: [`examples/311_step_gm_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/311_step_gm_agcu.jl)

### Phase Fractions vs Temperature — `StepResult`

Plot phase fractions by passing `Val(:phase_fraction)` as the second argument.

```julia
plot(result, Val(:phase_fraction))
```

This variant displays Phase 1 and Phase 2 fractions as separate series with `ylims=(0, 1)`.

| Attribute | Default |
|-----------|---------|
| xlabel | `"Temperature T [K]"` |
| ylabel | `"Phase Fraction"` |
| title | `"Phase Fractions vs Temperature"` |
| legend | `:topright` |
| ylims | `(0, 1)` |

See also: [`examples/313_step_twophase_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/313_step_twophase_agcu.jl)

### Binary Phase Diagram — `PhaseDiagramResult`

Plot a T-x phase diagram showing the miscibility gap boundary.

```julia
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")
result = map_phase_diagram(fcc, db, 800.0, 1300.0, 25.0)

plot(result)
```

The recipe filters converged points, constructs a closed boundary curve from left and right compositions, and fills the region between the two boundaries.

| Attribute | Default |
|-----------|---------|
| xlabel | `"Composition x"` |
| ylabel | `"Temperature T [K]"` |
| title | `"Binary Phase Diagram"` |
| fillalpha | 0.3 |

See also: [`examples/411_phase_diagram_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/411_phase_diagram_agcu.jl)

---

## Plots.jl Convenience Functions

These functions require `using Plots` and delegate to the corresponding recipes:

```julia
using OpenCALPHAD, Plots

plot_gibbs_curve(scan)         # equivalent to plot(scan::GridScanResult)
plot_phase_diagram(pd_result)  # equivalent to plot(result::PhaseDiagramResult)
```

See [`plot_phase_diagram`](@ref) and [`plot_gibbs_curve`](@ref) in the [API Reference](api.md).

## Publication-quality output (PythonPlot)

For publication-quality static figures (PDF / PNG via matplotlib), load the
PythonPlot extension. The API has **three layers** (L1 / L2 are unexported —
call them as `OpenCALPHAD.…`; L3 `savefig_publication` is exported):

| Layer | Function | Returns | Use |
|-------|----------|---------|-----|
| L3 | `savefig_publication(result, path; ...)` | `path` | save to PDF or PNG in one call (format from the file extension; delegates to L2) |
| L2 | `OpenCALPHAD.figure_publication(result; ...)` | `(fig, ax)` | a sized figure + axis to tweak before saving |
| L1 | `OpenCALPHAD.plot_on_axis!(ax, result; ...)` | `ax` | draw onto your own matplotlib axis (compose a subplot grid) |

```julia
using OpenCALPHAD, PythonPlot
```

On a headless machine or in CI (no display), select a non-interactive backend
before plotting: `PythonPlot.matplotlib.use("Agg")`.

In every layer the `result` argument may be a `GridScanResult` (Gibbs curve),
a `StepResult` (`kind = :gibbs` or `:phase_fraction`), or a `PhaseDiagramResult`
(binary phase diagram). The examples below use L3 (`savefig_publication`) for
each of these; the *Composing subplots* section further down shows L1 / L2.

### Gibbs energy curve — `GridScanResult`

```julia
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")
solver = GridSearchSolver(n_points = 51)
grid_result = scan_composition(fcc, 1000.0, db, solver)

savefig_publication(grid_result, "gibbs.pdf";
                    title = "Ag-Cu FCC Gibbs energy at T = 1000 K",
                    show_minimum = true)
```

`show_minimum = true` draws a vertical line at the composition of minimum G.

### `StepResult` — Gibbs vs axis or phase fractions

`kind` selects the plot. Default is `:gibbs`; `:phase_fraction` plots two
phase fractions with a forced `ylim = (0, 1)`.

```julia
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")
step_result = step_temperature(fcc, db, 0.5, 800.0, 1300.0, 25.0)

savefig_publication(step_result, "step_g.pdf")                       # kind = :gibbs (default)
savefig_publication(step_result, "step_pf.pdf"; kind = :phase_fraction)
```

### Binary phase diagram — `PhaseDiagramResult`

```julia
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")
pd_result = map_phase_diagram(fcc, db, 800.0, 1300.0, 25.0)

savefig_publication(pd_result, "phase_diagram.pdf";
                    axis_width_mm = 100.0, axis_height_mm = 70.0)
```

Points that did not converge are filtered out automatically; the region
between the two boundaries is shaded with `fill_color` at `fill_alpha`. If no
point converges, an empty axes with only the axis labels is drawn.

### Keyword reference

Keywords are documented on the functions below. The plot keywords for each
result type (and the label and title precedence) live on `plot_on_axis!`;
`figure_publication` adds the figure size keywords; `savefig_publication`
forwards everything:

```@docs
OpenCALPHAD.savefig_publication
OpenCALPHAD.figure_publication
OpenCALPHAD.plot_on_axis!
```

### Composing subplots with `plot_on_axis!` (L1) / `figure_publication` (L2)

Beyond `savefig_publication` (L3) shown above, the two lower layers let you
compose and tweak figures (see the layer table at the top of this section).

**L1 — compose multiple results into one figure** (e.g. a 2×2 poster panel).
You create the figure with matplotlib's `subplots`/`add_subplot` and own the
`close`:

```julia
using OpenCALPHAD, PythonPlot

fig = PythonPlot.figure(figsize = (10, 8))
OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 1), grid_result)
OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 2), step_result; kind = :gibbs)
OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 3), step_result; kind = :phase_fraction)
OpenCALPHAD.plot_on_axis!(fig.add_subplot(2, 2, 4), pd_result)
fig.savefig("composite.pdf")
PythonPlot.close(fig)            # caller owns the figure
```

The plot keyword arguments for each result type are the same as for `savefig_publication`
(e.g. `title`, `color`, `show_minimum`, `kind`); `plot_on_axis!` returns the
same `ax` for chaining.

**L2 — adjust before saving:**

```julia
fig, ax = OpenCALPHAD.figure_publication(grid_result; axis_width_mm = 100.0)
ax.axvline(0.5; linestyle = "--")     # add your own annotations
fig.savefig("annotated.pdf")
PythonPlot.close(fig)
```

> Note: inside the extension, figures are created with `PythonPlot.figure(...)`
> (not `subplots()`), but on the **caller** side `subplots()` is fine — you just
> own the resulting figure and must `close` it yourself.

### Installing PythonPlot

The extension is loaded automatically once both `OpenCALPHAD` and `PythonPlot`
are imported in the same session. Install PythonPlot with:

```julia
using Pkg
Pkg.add("PythonPlot")
```

CondaPkg.jl will install matplotlib on first use.

See also: [`examples/101_gibbs_energy_agcu_publication.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/101_gibbs_energy_agcu_publication.jl), [`examples/311_step_gm_agcu_publication.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/311_step_gm_agcu_publication.jl), [`examples/411_phase_diagram_agcu_publication.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/411_phase_diagram_agcu_publication.jl).

See [`savefig_publication`](@ref) in the [API Reference](api.md).

