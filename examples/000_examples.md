# OpenCALPHAD.jl Examples

Last-Modified: 2026-01-25

## Numbering System

| Range | Category | Description |
|-------|----------|-------------|
| 000 | Documentation | This index file |
| 1xx | Basic Calculations | Gibbs energy, chemical potential |
| 2xx | Equilibrium | Miscibility gap, two-phase equilibrium |
| 3xx | STEP / MAP | Temperature/composition sweeps |
| 4xx | Phase Diagrams | Complete phase diagrams |
| 5xx | Simulation Coupling | Driving force, diffusion potential |
| 6xx | Advanced | Magnetic contribution, parameter optimization |
| 9xx | Benchmarks | Validation against openCALPHAD (Fortran) |

**Material System Codes** (2nd digit):
| Code | System | TDB File |
|------|--------|----------|
| x1 | Ag-Cu | agcu.TDB |
| x2 | Al-Ni | alni-4slx.TDB (future) |

## File List

### 1xx: Basic Calculations

| File | Description |
|------|-------------|
| `101_gibbs_energy_agcu.jl` | Gibbs energy curve with common tangent construction |

### 2xx: Equilibrium

| File | Description |
|------|-------------|
| `211_miscibility_gap_agcu.jl` | Miscibility gap validation against literature |

### 3xx: STEP / Property Curves

| File | Description |
|------|-------------|
| `311_step_gm_agcu.jl` | G, H_mix, S_mix curves at fixed T (≈ step2.OCM) |
| `312_two_phase_crfe.jl` | Two-phase equilibrium: BCC gap + liquidus (≈ map2.OCM) |

### 4xx: Phase Diagrams

| File | Description |
|------|-------------|
| `411_phase_diagram_agcu.jl` | Ag-Cu binary phase diagram (≈ map1.OCM) |

### 6xx: Advanced

| File | Description |
|------|-------------|
| `611_magnetic_ihj.jl` | Inden-Hillert-Jarl magnetic contribution model |
| `621_parameter_optimization.jl` | Fitting R-K coefficients (+ combined figure for JOSS) |
| `621_parameter_optimization.ipynb` | Jupyter notebook version |
| `624_magnetic_optimization.jl` | Optimizing magnetic parameters (Tc, beta) |

### 9xx: Benchmarks / References

| File | Description |
|------|-------------|
| `901_validation_comparison.jl` | OpenCALPHAD.jl vs openCALPHAD (Fortran) comparison |
| `911_phase_diagram_validated.jl` | Phase diagram with Julia DSL validation |

### Notebooks

| File | Description |
|------|-------------|
| `411_phase_diagram_agcu.ipynb` | Interactive Jupyter notebook for Ag-Cu phase diagram |

## How to Run

```bash
# Navigate to the package directory
cd OpenCALPHAD.jl

# Run a single example
julia --project=. examples/101_gibbs_energy_agcu.jl

# Run with plots displayed
julia --project=. -i examples/411_phase_diagram_agcu.jl
```

## Output Images

Each example generates corresponding PNG files in the `examples/` directory:
- `101_gibbs_energy_agcu.png`
- `311_step_gm_agcu.png`
- `312_two_phase_crfe.png`
- `411_phase_diagram_agcu.png`
- `611_magnetic_ihj.png`
- `621_parameter_optimization.png`
- `621_parameter_optimization_combined.png` (JOSS paper figure)
- `624_magnetic_optimization.png`
- `901_validation_comparison.png`
- `911_phase_diagram_validated.png`

Note: `211_miscibility_gap_agcu.jl` is validation-only (no plot output).

## Dependencies

All examples use:
- `OpenCALPHAD` (this package)
- `Plots` (for visualization)

Some examples additionally require:
- `ForwardDiff` (for automatic differentiation)
- `Optim` (for parameter optimization)
- `JSON` (for reference data loading)

## See Also

- [Getting Started](../docs/src/getting_started.md) - Quick introduction
- [API Reference](../docs/src/api.md) - Function documentation
- [Julia DSL](../docs/src/julia_dsl.md) - Defining phases in Julia
