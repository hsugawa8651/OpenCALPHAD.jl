# OpenCALPHAD.jl

[![CI](https://github.com/hsugawa8651/OpenCALPHAD.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/hsugawa8651/OpenCALPHAD.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/hsugawa8651/OpenCALPHAD.jl/actions/workflows/Docs.yml/badge.svg)](https://hsugawa8651.github.io/OpenCALPHAD.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18670628.svg)](https://doi.org/10.5281/zenodo.18670628)

A pure Julia port of [OpenCALPHAD](https://github.com/sundmanbo/opencalphad) for CALPHAD (CALculation of PHAse Diagrams) thermodynamic calculations.

**Based on:** OpenCALPHAD (Fortran) commit [`e857325f`](https://github.com/sundmanbo/opencalphad/commit/e857325f5d38e32e41513186ddf1d8ae4b2e0349) (2026-01-01)

## Features (v0.1)

* Binary system equilibrium calculations (including miscibility gaps)
* TDB database parser
* Compound Energy Formalism (CEF) Gibbs energy calculation
* Phase diagram generation (STEP, MAP)
* Plots.jl integration

## Scope and Limitations (v0.1)

v0.1 supports **binary (2-component) systems** with the following phase structures:

| Structure | Example | Supported |
|-----------|---------|-----------|
| Single substitutional sublattice | `(A,B)₁` | ✅ |
| Substitutional + VA-only interstitial | `(A,B)₁(VA)ₙ` | ✅ |
| Substitutional + interstitial with VA | `(A,B)₁(C,VA)ₙ` | ✅ (VA fixed) |
| 3+ sublattices | `(A)₈(B)₄(C)₁₈` | ❌ |
| Multiple substitutional sublattices | `(A,B)₂(C)₁` | ❌ |
| No VA in interstitial sublattice | `(A,B)₁(C)ₙ` | ❌ |

**Key constraints:**
- First sublattice must be substitutional (contains the 2 components of interest)
- Second sublattice onwards must contain VA (vacancy), which is fixed to 1.0
- For multi-component databases (e.g., steel1.TDB), use `component_indices` parameter to select 2 components

### Available Functions

| Feature | Function | Description |
|---------|----------|-------------|
| Database I/O | `read_tdb()` | Read TDB thermodynamic database |
| Equilibrium | `calculate_equilibrium()` | Single-point equilibrium calculation |
| Composition scan | `scan_composition()` | G vs x calculation (supports `component_indices`) |
| Temperature sweep | `step_temperature()` | STEP calculation (supports `component_indices`) |
| Miscibility gap | `find_miscibility_gap()` | Find two-phase region (supports `component_indices`) |
| Phase diagram | `map_phase_diagram()` | MAP calculation (supports `component_indices`) |
| Plotting | `plot()` | Plot results (requires Plots.jl) |

## Installation

```julia
using Pkg
Pkg.add("OpenCALPHAD")
```

## Quick Example

```julia
using OpenCALPHAD
using Plots

# Load thermodynamic database
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")

# Calculate binary phase diagram
result = map_phase_diagram(fcc, db, 600.0, 1200.0, 25.0)

# Plot
plot(result)
```

## Documentation

Documentation is available [here](https://hsugawa8651.github.io/OpenCALPHAD.jl/dev). Build locally:

```bash
cd docs
julia --project make.jl
```

## Validation

OpenCALPHAD.jl is validated against openCALPHAD (Fortran) through 175 reference tests.

```bash
# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Run reference tests separately
cd reftest
julia --project run_reftest.jl
```

See `reftest/README.md` for details.

## Citation

If you use OpenCALPHAD.jl in your research, please cite both:

**OpenCALPHAD.jl:**
> Sugawara, H. (2026). OpenCALPHAD.jl: Julia port of openCALPHAD for CALPHAD thermodynamic calculations (Version 0.1.0) [Computer software]. https://doi.org/10.5281/zenodo.18670628

**Original openCALPHAD:**
> Sundman, B., Kattner, U. R., Palumbo, M., & Fries, S. G. (2015). OpenCalphad - a free thermodynamic software. *Integrating Materials and Manufacturing Innovation*, 4:1. [doi:10.1186/s40192-014-0029-1](https://doi.org/10.1186/s40192-014-0029-1)

## Contributing

Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/hsugawa8651/OpenCALPHAD.jl/issues).
Before opening a pull request, start an issue or a discussion on the topic.
This project follows the [Julia Community Standards](https://julialang.org/community/standards/).

## References

* OpenCALPHAD (Fortran): [github.com/sundmanbo/opencalphad](https://github.com/sundmanbo/opencalphad)
* Sundman et al., "OpenCalphad - a free thermodynamic software", *Integrating Materials and Manufacturing Innovation* (2015) 4:1

## License

GPL-3.0-or-later (same as OpenCALPHAD Fortran)

See [LICENSE](LICENSE) for details.
