# Julia DSL

In addition to reading TDB files, OpenCALPHAD.jl provides a Julia-native Domain Specific Language (DSL)
for defining thermodynamic models directly in Julia code.

!!! info "Advantage of Julia DSL"
    The Julia DSL enables **automatic differentiation** through iterative solvers,
    making gradient-based parameter optimization possible.

!!! tip "Example"
    See [`examples/131_julia_dsl_agcu.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/131_julia_dsl_agcu.jl) for a complete working example.

## Basic Usage

### Creating a Phase

```julia
using OpenCALPHAD

# Create a binary FCC phase with single sublattice
fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])
#            name     sites   constituents
```

- `name`: Phase identifier
- `sites`: Sublattice site multiplicities (e.g., `[1.0]` for single sublattice)
- `constituents`: Species on each sublattice

### Setting Gibbs Energy (`set_G!`)

Define the Gibbs energy of pure components (endmembers):

```julia
# GHSER for pure Ag (from SGTE database)
function ghser_ag(T)
    return -7209.512 + 118.200733*T - 23.8463314*T*log(T) -
           0.001790585*T^2 - 3.98587e-7*T^3 - 12011/T
end

# GHSER for pure Cu
function ghser_cu(T)
    return -7770.458 + 130.485403*T - 24.112392*T*log(T) -
           0.00265684*T^2 + 1.29223e-7*T^3 + 52478/T
end

# Set endmember Gibbs energies
set_G!(fcc, [:AG], ghser_ag)
set_G!(fcc, [:CU], ghser_cu)
```

### Setting Interaction Parameters (`set_L!`)

Define Redlich-Kister interaction parameters:

```julia
# L₀ (0th order interaction)
set_L!(fcc, [:AG, :CU], 0, T -> 33819.09 - 8.1236*T)

# L₁ (1st order interaction)
set_L!(fcc, [:AG, :CU], 1, T -> -5601.87 + 1.32997*T)

# L₂ (2nd order, if needed)
# set_L!(fcc, [:AG, :CU], 2, T -> 1000.0)
```

The excess Gibbs energy is computed as:

$$G^{ex} = x_{Ag} x_{Cu} \sum_{\nu=0}^{n} L_{\nu} (x_{Ag} - x_{Cu})^{\nu}$$

## Comparison: TDB vs Julia DSL

TDB (Thermo-Calc Database) is a widely-used text format for storing thermodynamic parameters.
It originated from [Thermo-Calc](https://thermocalc.com/) and is supported by most CALPHAD software
including [openCALPHAD](https://github.com/sundmanbo/opencalphad).

### TDB Format (openCALPHAD)

```
FUNCTION GHSERAG 298.15 -7209.512+118.200733*T-23.8463314*T*LN(T)
    -.001790585*T**2-3.98587E-07*T**3-12011*T**(-1); 1235.08 Y
    -15095.314+190.265169*T-33.472*T*LN(T)+1.412186E+29*T**(-9); 3000 N !

PARAMETER G(FCC_A1,AG:VA;0) 298.15 +GHSERAG#; 3000 N !
PARAMETER G(FCC_A1,CU:VA;0) 298.15 +GHSERCU#; 3200 N !
PARAMETER G(FCC_A1,AG,CU:VA;0) 298.15 +33819.1-8.1236*T; 6000 N !
PARAMETER G(FCC_A1,AG,CU:VA;1) 298.15 -5601.9+1.32997*T; 6000 N !
```

### Julia DSL (OpenCALPHAD.jl)

```julia
fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])

set_G!(fcc, [:AG], T -> -7209.512 + 118.200733*T - 23.8463314*T*log(T) -
                        0.001790585*T^2 - 3.98587e-7*T^3 - 12011/T)
set_G!(fcc, [:CU], T -> -7770.458 + 130.485403*T - 24.112392*T*log(T) -
                        0.00265684*T^2 + 1.29223e-7*T^3 + 52478/T)
set_L!(fcc, [:AG, :CU], 0, T -> 33819.09 - 8.1236*T)
set_L!(fcc, [:AG, :CU], 1, T -> -5601.87 + 1.32997*T)
```

| Feature | TDB | Julia DSL |
|---------|-----|-----------|
| Syntax | CALPHAD-specific | Native Julia |
| Temperature ranges | `Y`/`N` markers | `if/else` |
| Custom functions | Limited | Arbitrary Julia code |
| Auto-differentiation | No | Yes (ForwardDiff) |
| External libraries | No | Yes |

## Advanced Features

### Conditional Logic (Temperature Ranges)

TDB uses `Y`/`N` markers for temperature ranges. Julia DSL uses native `if/else`:

```julia
function ghser_ag(T)
    if T < 1235.08
        return -7209.512 + 118.200733*T - 23.8463314*T*log(T) -
               0.001790585*T^2 - 3.98587e-7*T^3 - 12011/T
    else
        return -15095.314 + 190.265169*T - 33.472*T*log(T) +
               1.412186e29*T^(-9)
    end
end
```

### Magnetic Contribution (Inden-Hillert-Jarl)

Complex models with conditional branching are naturally expressed.

Parameters:
- `Tc`: Curie temperature [K]
- `beta` (β): Average magnetic moment per atom [Bohr magneton]
- `p`: Structure-dependent coefficient (BCC: 0.4, FCC: 0.28)

```julia
function G_magnetic(T, Tc, beta, p)
    tau = T / Tc
    A = 518/1125 + (11692/15975) * (1/p - 1)

    if tau < 1  # Ferromagnetic
        g = 1 - (79/(140*p)*tau^(-1) + ...) / A
    else        # Paramagnetic
        g = -(tau^(-5)/10 + tau^(-15)/315 + tau^(-25)/1500) / A
    end

    return R * T * log(beta + 1) * g
end

# Use in phase definition
set_G!(bcc_fe, [:FE], T -> ghser_fe(T) + G_magnetic(T, 1043.0, 2.22, 0.4))
```

See [`examples/611_magnetic_ihj.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/611_magnetic_ihj.jl)

### External Library Integration

Julia DSL allows calling external libraries:

```julia
using Interpolations

# Load experimental data
T_data = [300, 400, 500, 600, 700]
G_data = [-1000, -2000, -3500, -5000, -6500]
itp = linear_interpolation(T_data, G_data)

# Use in Gibbs energy
set_G!(phase, [:X], T -> itp(T))
```

## Calculation without Database

With Julia DSL, no `Database` argument is needed:

```julia
# TDB workflow (requires db)
gap = find_miscibility_gap(fcc, 1000.0, db)

# Julia DSL workflow (no db)
gap = find_miscibility_gap(fcc, 1000.0)
```

## Parameter Optimization

The key advantage of Julia DSL is enabling automatic differentiation:

```julia
using ForwardDiff
using Optim

function loss(params)
    L0, L1 = params
    fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])
    set_G!(fcc, [:AG], ghser_ag)
    set_G!(fcc, [:CU], ghser_cu)
    set_L!(fcc, [:AG, :CU], 0, T -> L0[1] + L0[2]*T)
    set_L!(fcc, [:AG, :CU], 1, T -> L1[1] + L1[2]*T)

    error = 0.0
    for (T, x_exp) in experimental_data
        gap = find_miscibility_gap(fcc, T)
        error += (gap.x1 - x_exp)^2
    end
    return error
end

# Gradient via ForwardDiff (through Newton iterations!)
grad = ForwardDiff.gradient(loss, initial_params)

# Optimize with BFGS
result = optimize(loss, loss_grad!, initial_params, BFGS())
```

See [`examples/621_parameter_optimization.jl`](https://github.com/hsugawa8651/OpenCALPHAD.jl/blob/main/examples/621_parameter_optimization.jl)

## API Reference

| Function | Description |
|----------|-------------|
| `Phase(name, sites, constituents)` | Create a new phase |
| `set_G!(phase, species, func)` | Set endmember Gibbs energy |
| `set_L!(phase, species, order, func)` | Set interaction parameter |
| `has_julia_parameters(phase)` | Check if phase has Julia DSL parameters |
| `validate_ad_compatibility(func)` | Verify function works with ForwardDiff |

## See Also

- [Workflow](workflow.md) - Overall calculation workflow
- [Getting Started](getting_started.md) - Installation and basic usage
- [API Reference](api.md) - Full function documentation
