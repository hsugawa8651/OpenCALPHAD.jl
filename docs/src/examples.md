# Examples

This page provides an overview of all examples included with OpenCALPHAD.jl.

## Example Categories

| Range | Category | Description |
|-------|----------|-------------|
| 1xx | Basic Calculations | Gibbs energy, chemical potential |
| 2xx | Equilibrium | Miscibility gap, two-phase equilibrium |
| 3xx | STEP Calculations | Temperature stepping with phase detection |
| 4xx | Phase Diagrams | Binary phase diagrams |
| 6xx | Advanced | Magnetic models, parameter optimization |
| 9xx | Validation | Comparison with openCALPHAD (Fortran) |

## Example List

### Basic Calculations (1xx)

| File | Description |
|------|-------------|
| [`101_gibbs_energy_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/101_gibbs_energy_agcu.jl) | Gibbs energy curve with common tangent construction |
| [`102_gibbs_energy_crmo.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/102_gibbs_energy_crmo.jl) | Cr-Mo from multi-component database (manual y-matrix) |
| [`122_scan_crmo.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/122_scan_crmo.jl) | Cr-Mo with `component_indices` (high-level API) |
| [`131_julia_dsl_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/131_julia_dsl_agcu.jl) | **Julia DSL basics** - define phases without TDB |
| [`131_julia_dsl_agcu.ipynb`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/131_julia_dsl_agcu.ipynb) | Jupyter notebook version |

### Equilibrium (2xx)

| File | Description |
|------|-------------|
| [`211_miscibility_gap_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/211_miscibility_gap_agcu.jl) | Miscibility gap validation against literature |

### STEP Calculations (3xx)

| File | Description |
|------|-------------|
| [`311_step_gm_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/311_step_gm_agcu.jl) | Ag-Cu STEP calculation with Gm plot |
| [`311_step_gm_agcu.ipynb`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/311_step_gm_agcu.ipynb) | Jupyter notebook version (interactive) |
| [`312_two_phase_crfe.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/312_two_phase_crfe.jl) | Cr-Fe two-phase equilibrium |
| [`313_step_twophase_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/313_step_twophase_agcu.jl) | Two-phase detection with step_temperature |
| [`313_step_twophase_agcu.ipynb`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/313_step_twophase_agcu.ipynb) | Jupyter notebook version (interactive) |
| [`322_step_crmo.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/322_step_crmo.jl) | Cr-Mo with `component_indices` |

### Phase Diagrams (4xx)

| File | Description |
|------|-------------|
| [`411_phase_diagram_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/411_phase_diagram_agcu.jl) | Ag-Cu binary phase diagram |
| [`411_phase_diagram_agcu.ipynb`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/411_phase_diagram_agcu.ipynb) | Jupyter notebook version |

### Advanced (6xx)

| File | Description |
|------|-------------|
| [`611_magnetic_ihj.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/611_magnetic_ihj.jl) | Inden-Hillert-Jarl magnetic contribution model |
| [`621_parameter_optimization.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/621_parameter_optimization.jl) | Fitting Redlich-Kister coefficients |
| [`621_parameter_optimization.ipynb`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/621_parameter_optimization.ipynb) | Jupyter notebook version |
| [`624_magnetic_optimization.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/624_magnetic_optimization.jl) | Optimizing magnetic parameters (Tc, β) |

### Validation (9xx)

| File | Description |
|------|-------------|
| [`901_validation_comparison.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/901_validation_comparison.jl) | OpenCALPHAD.jl vs openCALPHAD (Fortran) |
| [`911_phase_diagram_validated.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/911_phase_diagram_validated.jl) | Phase diagram with Julia DSL validation |

## Running Examples

```bash
# Navigate to the package directory
cd OpenCALPHAD.jl

# Run a single example
julia --project=. examples/101_gibbs_energy_agcu.jl

# Run interactively (keeps plot window open)
julia --project=. -i examples/411_phase_diagram_agcu.jl
```

## See Also

- [Getting Started](getting_started.md) - Quick introduction with code snippets
- [Julia DSL](julia_dsl.md) - Defining phases directly in Julia
- [API Reference](api.md) - Complete function documentation
