# Plotting

OpenCALPHAD.jl provides three plotting paths via package extensions:

| Extension | Trigger | Functionality |
|-----------|---------|---------------|
| RecipesBaseExt | `using RecipesBase` (or any backend) | `plot(result)` recipes for all result types |
| PlotsExt | `using Plots` | Convenience functions `plot_gibbs_curve`, `plot_phase_diagram` |
| PythonPlotExt | `using PythonPlot` | Planned |

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
result = step_calculation(fcc, db, 300.0, 1400.0, 10.0; x_overall=0.5)

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

The recipe filters converged points, constructs a closed boundary curve from left and right compositions, and fills the two-phase region.

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

