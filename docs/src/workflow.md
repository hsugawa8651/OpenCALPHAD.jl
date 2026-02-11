# Workflow

## Overview

OpenCALPHAD.jl provides two main workflows for thermodynamic calculations:

1. **TDB Workflow**: Load parameters from TDB database files
2. **Julia DSL Workflow**: Define parameters directly in Julia code

```
┌─────────────────────────────────────────────────────────────────┐
│                      TDB Workflow                               │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   agcu.TDB ──► read_tdb() ──► Database                         │
│                                   │                             │
│                              get_phase()                        │
│                                   │                             │
│                                   ▼                             │
│                                 Phase ──────────────────────┐   │
│                                                             │   │
└─────────────────────────────────────────────────────────────│───┘
                                                              │
┌─────────────────────────────────────────────────────────────│───┐
│                    Julia DSL Workflow                       │   │
├─────────────────────────────────────────────────────────────│───┤
│                                                             │   │
│   Phase("FCC_A1", ...) ──► set_G!() ──► set_L!() ──────────┘   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                     Equilibrium Calculation                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   find_miscibility_gap()  ──► MiscibilityGapResult             │
│   map_phase_diagram()     ──► PhaseDiagramResult               │
│   step_temperature()      ──► StepResult                       │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Analysis / Output                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   plot()              ──► Phase diagram visualization           │
│   ForwardDiff         ──► Parameter optimization               │
│   chemical_potential  ──► Phase-field coupling                  │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Basic Workflow (TDB)

### Step-by-Step

1. Load TDB database with `read_tdb()`
2. Get phase with `get_phase()`
3. Calculate equilibrium with `find_miscibility_gap()` or `calculate_equilibrium()`
4. Plot results with `plot()`

### Example: Ag-Cu Phase Diagram

See also: [`examples/plot_phase_diagram.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/plot_phase_diagram.jl)

```julia
using OpenCALPHAD
using Plots

# 1. Load database
db = read_tdb("agcu.TDB")

# 2. Get phase
fcc = get_phase(db, "FCC_A1")

# 3. Calculate phase diagram
result = map_phase_diagram(fcc, db, 800.0, 1300.0, 25.0)

# 4. Plot
plot(result)
```

### Example: Single Temperature Equilibrium

```julia
# Find miscibility gap at 1000 K
gap = find_miscibility_gap(fcc, 1000.0, db)

println("Ag-rich boundary: x(Ag) = $(gap.x1)")
println("Cu-rich boundary: x(Ag) = $(gap.x2)")
```

## Julia DSL Workflow

The Julia DSL allows defining thermodynamic parameters directly as Julia functions,
enabling automatic differentiation for parameter optimization.

### Step-by-Step

1. Create `Phase` with constituents
2. Set Gibbs energy with `set_G!()`
3. Set interaction parameters with `set_L!()`
4. Calculate equilibrium (no database argument needed)

### Example: Ag-Cu with Julia DSL

See also: [`examples/plot_phase_diagram_validated.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/plot_phase_diagram_validated.jl)

```julia
using OpenCALPHAD

# GHSER functions (from SGTE database)
function ghser_ag(T)
    return -7209.512 + 118.200733*T - 23.8463314*T*log(T) -
           0.001790585*T^2 - 3.98587e-7*T^3 - 12011/T
end

function ghser_cu(T)
    return -7770.458 + 130.485403*T - 24.112392*T*log(T) -
           0.00265684*T^2 + 1.29223e-7*T^3 + 52478/T
end

# 1. Create phase
fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])

# 2. Set Gibbs energies
set_G!(fcc, [:AG], ghser_ag)
set_G!(fcc, [:CU], ghser_cu)

# 3. Set interaction parameters (Redlich-Kister)
set_L!(fcc, [:AG, :CU], 0, T -> 33819.09 - 8.1236*T)   # L₀
set_L!(fcc, [:AG, :CU], 1, T -> -5601.87 + 1.32997*T)  # L₁

# 4. Calculate (no db argument)
gap = find_miscibility_gap(fcc, 1000.0)
```

## Parameter Optimization Workflow

The Julia DSL enables gradient-based optimization of thermodynamic parameters
using automatic differentiation through iterative solvers.

### Step-by-Step

1. Define phase with `set_G!()` and `set_L!()`
2. Define loss function comparing calculated vs. experimental data
3. Optimize with `Optim.jl` (BFGS or other gradient-based method)

### Example: Fitting Redlich-Kister Parameters

See also: [`examples/parameter_optimization_combined.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/parameter_optimization_combined.jl)

```julia
using OpenCALPHAD
using Optim
using ForwardDiff

# Experimental data: (T, x_low, x_high)
experimental_data = [
    (800.0,  0.038, 0.991),
    (1000.0, 0.103, 0.966),
    (1200.0, 0.224, 0.905),
]

# Create phase with parameterized L
function create_phase(L0_a, L0_b, L1_a, L1_b)
    fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])
    set_G!(fcc, [:AG], ghser_ag)
    set_G!(fcc, [:CU], ghser_cu)
    set_L!(fcc, [:AG, :CU], 0, T -> L0_a + L0_b * T)
    set_L!(fcc, [:AG, :CU], 1, T -> L1_a + L1_b * T)
    return fcc
end

# Loss function
function loss(params)
    L0_a, L0_b, L1_a, L1_b = params
    fcc = create_phase(L0_a, L0_b, L1_a, L1_b)
    error = 0.0
    for (T, x_low, x_high) in experimental_data
        gap = find_miscibility_gap(fcc, T)
        if !isnothing(gap)
            error += (gap.x1 - x_low)^2 + (gap.x2 - x_high)^2
        end
    end
    return error
end

# Gradient function (ForwardDiff)
function loss_grad!(G, params)
    G .= ForwardDiff.gradient(loss, params)
end

# Optimize
initial_params = [35000.0, -10.0, -5000.0, 1.0]
result = optimize(loss, loss_grad!, initial_params, BFGS())

println("Optimized parameters: ", Optim.minimizer(result))
```

## Magnetic Contribution (Advanced)

The Inden-Hillert-Jarl model for magnetic contributions uses conditional branching
(ferromagnetic vs. paramagnetic regions). ForwardDiff handles this automatically.

See also: [`examples/magnetic_contribution.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/magnetic_contribution.jl)

```julia
function G_magnetic(T, Tc, beta, p)
    tau = T / Tc
    if tau < 1  # Ferromagnetic
        g = 1 - (79/(140*p)*tau^(-1) + ...) / A
    else        # Paramagnetic
        g = -(tau^(-5)/10 + ...) / A
    end
    return R * T * log(beta + 1) * g
end

# Use in DSL
set_G!(bcc_fe, [:FE], T -> ghser_fe(T) + G_magnetic(T, 1043.0, 2.22, 0.4))
```

## Summary

| Workflow | Use Case | Key Functions |
|----------|----------|---------------|
| TDB | Standard calculations | `read_tdb()`, `get_phase()` |
| Julia DSL | Custom models, optimization | `Phase()`, `set_G!()`, `set_L!()` |
| Optimization | Parameter fitting | `ForwardDiff.gradient()`, `optimize()` |

## See Also

- [Getting Started](getting_started.md) - Installation and basic usage
- [API Reference](api.md) - Full function documentation
