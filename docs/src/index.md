# OpenCALPHAD.jl

A pure Julia implementation of CALPHAD (CALculation of PHAse Diagrams) thermodynamic calculations.

## Overview

OpenCALPHAD.jl provides tools for:

- Reading thermodynamic databases (TDB format)
- Calculating Gibbs energy using the Compound Energy Formalism (CEF)
- Finding phase equilibria including miscibility gaps
- Generating binary phase diagrams
- Plotting results with Plots.jl integration

## Features

- **Pure Julia**: No external dependencies on Fortran or C libraries
- **Automatic Differentiation**: Uses ForwardDiff.jl for accurate gradient calculations
- **Julia DSL**: Define thermodynamic models directly as Julia functions
- **Phase Field Integration**: Utility functions for coupling with phase field simulations
- **Modern API**: Designed for the Julia ecosystem
- **Validated**: Reference tests against openCALPHAD (Fortran)

## Scope and Limitations (v0.1)

v0.1 is designed for **binary (2-component) systems** with specific phase structures.

### Supported Phase Structures

| Structure | Description | Example |
|-----------|-------------|---------|
| `(A,B)₁` | Single substitutional sublattice | LIQUID |
| `(A,B)₁(VA)ₙ` | Substitutional + VA-only interstitial | FCC_A1 in Ag-Cu |
| `(A,B,...)₁(C,VA)ₙ` | Substitutional + interstitial with VA | BCC_A2 in steel1.TDB |

### Not Supported in v0.1

| Structure | Reason | Planned |
|-----------|--------|---------|
| 3+ sublattices (e.g., SIGMA) | Multiple substitutional sites | v0.2 |
| 3+ component systems | Requires composition vector | v0.3 |
| Interstitial without VA | Cannot fix interstitial site | Future |

### Key Constraints

1. **First sublattice**: Must be substitutional and contain the 2 components of interest
2. **Second sublattice onwards**: Must contain VA (vacancy), fixed to y(VA)=1.0
3. **Multi-component databases**: Use `component_indices` parameter to select which 2 components to calculate

### Example: Cr-Mo from steel1.TDB

```julia
db = read_tdb("steel1.TDB")  # 6-component database
bcc = get_phase(db, "BCC_A2")
# BCC_A2: (CR,FE,MO,SI,V)₁(C,VA)₃
# CR=index 1, MO=index 3

# Select Cr-Mo binary from 6-component database
result = step_temperature(bcc, db, 0.5, 1500.0, 2500.0, 50.0;
                          component_indices=(1, 3))  # (CR, MO)
```

## Installation

```julia
using Pkg
Pkg.add("OpenCALPHAD")
```

Or in the Pkg REPL:

```
]add OpenCALPHAD
```

## Quick Example

```julia
using OpenCALPHAD
using Plots

# Load thermodynamic database
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")

# Calculate phase diagram
result = map_phase_diagram(fcc, db, 600.0, 1200.0, 25.0)

# Plot
plot(result)
```

## References

- OpenCALPHAD (Fortran): [github.com/sundmanbo/opencalphad](https://github.com/sundmanbo/opencalphad)
- Sundman et al., "OpenCalphad - a free thermodynamic software", Integrating Materials and Manufacturing Innovation (2015) 4:1

## License

GPL-3.0-or-later (same as OpenCALPHAD Fortran)

This software is a derivative work of openCALPHAD and is distributed under the
[GNU General Public License v3.0 or later](https://www.gnu.org/licenses/gpl-3.0.html).
