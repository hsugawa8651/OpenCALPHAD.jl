# Getting Started

This guide covers the basic usage of OpenCALPHAD.jl for thermodynamic calculations.

## Reading a Database

Thermodynamic data is read from TDB (Thermo-Calc Database) files:

```julia
using OpenCALPHAD

db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")
```

## Calculating Gibbs Energy

Calculate the Gibbs energy at a specific temperature and composition:

```julia
y = [0.7 0.3]  # Site fractions: 70% Ag, 30% Cu
T = 1000.0    # Temperature [K]
G = calculate_gibbs_energy(fcc, T, y, db)
```

See also: [`101_gibbs_energy_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/101_gibbs_energy_agcu.jl)

## Finding Miscibility Gap

For phases with miscibility gaps (like FCC in Ag-Cu):

```julia
gap = find_miscibility_gap(fcc, 1000.0, db)
# gap.x1 = Cu-rich boundary, gap.x2 = Ag-rich boundary
```

See also: [`211_miscibility_gap_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/211_miscibility_gap_agcu.jl)

## Phase Diagram Calculation

Calculate phase boundaries over a temperature range:

```julia
result = map_phase_diagram(fcc, db, 600.0, 1200.0, 25.0)

# Access boundary compositions
T = temperatures(result)
x_left = boundaries(result, 1)   # Cu-rich boundary
x_right = boundaries(result, 2)  # Ag-rich boundary
```

See also: [`411_phase_diagram_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/411_phase_diagram_agcu.jl)

## Plotting

OpenCALPHAD.jl integrates with Plots.jl through recipe functions:

```julia
using Plots

result = map_phase_diagram(fcc, db, 600.0, 1200.0, 25.0)
plot(result)
savefig("phase_diagram.png")
```

## State Variables

| Symbol | Description | Unit |
|--------|-------------|------|
| T | Temperature | K |
| P | Pressure | Pa |
| G | Gibbs energy | J/mol |
| x | Mole fraction | - |
| y | Site fraction | - |

## Quick Example

```julia
using OpenCALPHAD
using Plots

# Load database and get phase
db = read_tdb("agcu.TDB")
fcc = get_phase(db, "FCC_A1")

# Calculate and plot phase diagram
result = map_phase_diagram(fcc, db, 600.0, 1200.0, 25.0)
plot(result, title="Ag-Cu Phase Diagram")
savefig("agcu_phase_diagram.png")
```

## Next Steps

- [Examples](examples.md) - Complete list of runnable examples
- [Julia DSL](julia_dsl.md) - Define thermodynamic models in Julia
- [API Reference](api.md) - Complete function documentation
