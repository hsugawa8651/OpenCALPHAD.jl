# API Reference

## Database

```@docs
read_tdb
Database
get_phase
get_element
get_function
```

## Types

```@docs
Phase
Element
Parameter
JuliaParameter
GFunction
Equilibrium
PhaseResult
ComponentResult
EquilCondition
ConditionSet
ConditionType
```

## Gibbs Energy

```@docs
calculate_gibbs_energy
```

## Thermodynamic Properties

Functions for calculating thermodynamic properties (entropy, enthalpy, mixing properties).
Currently only supports binary (2-component) systems.

```@docs
calculate_entropy
calculate_enthalpy
calculate_mixing_gibbs
calculate_mixing_enthalpy
calculate_mixing_entropy
thermodynamic_properties
```

## Solvers

```@docs
DiffMethod
EquilibriumSolver
GridSearchSolver
NewtonSolver
HybridSolver
SinglePhaseResult
GridScanResult
scan_composition
scan_composition_multi
minimize_gibbs
find_minimum_envelope
find_two_phase_region
clamp_site_fractions!
```

## Two-Phase Equilibrium

```@docs
find_miscibility_gap
calculate_miscibility_equilibrium
find_common_tangent
calculate_two_phase_equilibrium
find_phase_boundaries
MiscibilityGapResult
TwoPhaseResult
```

## STEP Calculation

```@docs
step_temperature
StepResult
StepPoint
temperatures
gibbs_energies
phase_fractions
```

## MAP Calculation

```@docs
map_phase_diagram
PhaseDiagramResult
PhaseDiagramPoint
boundaries
```

## Differentiation

```@docs
ForwardDiffMethod
FiniteDiffMethod
compute_gradient
compute_hessian
compute_chemical_potential
```

## Julia DSL

Functions for defining thermodynamic models directly in Julia code.

```@docs
set_G!
set_L!
has_julia_parameters
validate_ad_compatibility
is_endmember
is_interaction
get_julia_endmembers
get_julia_interactions
```

### Database Building (DSL)

```@docs
add_element!
add_elements!
add_phase!
```

## Phase Field Integration

Utility functions for coupling with phase field simulations.

```@docs
driving_force
chemical_potential
diffusion_potential
phase_field_params
```

## Equilibrium Calculation

```@docs
calculate_equilibrium
step_calculation
set_condition!
degrees_of_freedom
stable_phases
n_stable_phases
```

## Utility Functions

```@docs
num_sublattices
total_sites
evaluate_gfunction
compile_gfunction!
evaluate_parameter
```

## Plotting

```@docs
plot_phase_diagram
plot_gibbs_curve
```

## Index

```@index
```
