---
title: 'OpenCALPHAD.jl: A Pure Julia Implementation of CALPHAD Thermodynamic Calculations'
tags:
  - Julia
  - thermodynamics
  - CALPHAD
  - phase diagrams
  - materials science
authors:
  - name: Hiroharu Sugawara
    orcid: 0000-0002-0071-2396
    corresponding: true
    email: hsugawa@tmu.ac.jp
    affiliation: 1
affiliations:
  - name: Graduate School of Systems Design, Tokyo Metropolitan University, Japan
    index: 1
date: 23 January 2026
bibliography: paper.bib
---

# Summary

Phase diagrams map the thermodynamically stable phases of materials as functions of temperature, pressure, and composition. They are essential tools in materials science and engineering, providing the foundation for alloy design, heat treatment optimization, and predicting material behavior under different processing conditions. The CALculation of PHAse Diagrams (CALPHAD) method enables computational prediction of these diagrams from thermodynamic databases, supporting the development of new materials for applications ranging from aerospace alloys to battery electrodes.

OpenCALPHAD.jl is a pure Julia port of openCALPHAD [@Sundman2015], an open-source CALPHAD software developed by Bo Sundman, a co-developer of Thermo-Calc and one of the pioneers of computational thermodynamics. The package provides tools for reading thermodynamic databases in TDB format (a widely-used text format originating from Thermo-Calc), calculating Gibbs energies using the Compound Energy Formalism—a standard model for describing the energetics of solid solutions with multiple sublattices—and computing phase equilibria including miscibility gaps in binary systems.

Unlike existing CALPHAD software that relies on Fortran or C backends, OpenCALPHAD.jl is written entirely in Julia, eliminating the need for external binary dependencies and simplifying installation across platforms. The package leverages Julia's automatic differentiation capabilities through ForwardDiff.jl [@Revels2016], enabling accurate and efficient computation of thermodynamic derivatives required for equilibrium calculations.

OpenCALPHAD.jl integrates seamlessly with the Julia ecosystem, including Plots.jl for visualization of phase diagrams and Gibbs energy curves. The implementation has been validated against openCALPHAD (Fortran) v6.100 [@Sundman2015] for binary Ag-Cu system calculations, demonstrating numerical agreement within 0.02% for phase boundaries.

The package is designed with educational use in mind, featuring readable code that closely follows the mathematical formulation of CALPHAD. This makes it suitable for teaching computational thermodynamics and for researchers who wish to understand or extend the underlying algorithms.

# Statement of Need

The CALPHAD method [@Lukas2007] is the standard approach for calculating phase diagrams and thermodynamic properties in materials science. However, existing implementations present barriers to adoption:

- **Commercial tools** (Thermo-Calc, FactSage) require expensive licenses
- **Open-source Fortran codes** (openCALPHAD) require Fortran compilers and manual implementation of derivatives
- **Python implementations** (pycalphad [@Otis2017]) depend on external C++ libraries (SymEngine) for symbolic differentiation

| | Thermo-Calc | openCALPHAD | pycalphad | OpenCALPHAD.jl |
|---|---|---|---|---|
| License | Commercial | GPL-3.0 | MIT | GPL-3.0 |
| Language | Proprietary | Fortran | Python + C++ | Pure Julia |
| Differentiation | Analytical | Analytical | Symbolic | AD |
| External dependencies | — | Fortran compiler | SymEngine (C++) | None |

OpenCALPHAD.jl addresses these limitations:

1. **Free and open-source**: Licensed under GPL-3.0, freely available for academic and commercial use.

2. **Pure Julia implementation** [@Bezanson2017]: No compilation required, no external binary dependencies. All code—from TDB parsing to equilibrium solving—is written in Julia, enabling users to inspect and understand any component.

3. **Automatic differentiation**: ForwardDiff.jl [@Revels2016] computes exact derivatives through arbitrary Julia code, including iterative solvers—unlike symbolic differentiation, which cannot differentiate through iterative algorithms.

Target users include materials scientists performing thermodynamic calculations, educators teaching computational thermodynamics, and developers building thermodynamic modeling workflows in Julia.

# Features

## Pure Julia Implementation

In addition to reading TDB files, OpenCALPHAD.jl provides a Julia DSL that allows users to define thermodynamic models directly as Julia functions:

```julia
fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])
set_G!(fcc, [:AG], T -> -7209 + 118*T - 24*T*log(T))  # Gibbs energy of pure Ag
set_G!(fcc, [:CU], T -> -7770 + 130*T - 24*T*log(T))  # Gibbs energy of pure Cu
set_L!(fcc, [:AG, :CU], 0, T -> 35620 - 1.85*T)       # Interaction parameter

gap = find_miscibility_gap(fcc, 1000.0)  # Returns phase boundaries at 1000 K
```

## Automatic Differentiation

Thermodynamic parameters can be fitted to experimental phase boundary data. Because ForwardDiff.jl can differentiate through iterative equilibrium solvers, this fitting can use gradient-based optimization methods like BFGS (\autoref{fig:optimization}).

![Fitting Redlich-Kister interaction parameters ($L_0$, $L_1$) to target phase boundary data. (a) Miscibility gap evolution from initial guess (red) to final convergence (blue). Black stars show target boundaries. (b) Sum of squared errors showing 5 orders of magnitude improvement.\label{fig:optimization}](figures/parameter_optimization_combined.png){ width=100% }

## Validation

The package has been validated for binary phase diagram calculations against the reference Fortran implementation openCALPHAD v6.100 [@Sundman2015], demonstrating numerical agreement within 0.02% for phase boundaries (\autoref{fig:phase_diagram}).

![Ag-Cu FCC miscibility gap calculated using OpenCALPHAD.jl (lines) with thermodynamic parameters from the COST 531 database. Red circles show reference data from openCALPHAD (Fortran) v6.100, confirming numerical agreement within 0.02%.\label{fig:phase_diagram}](figures/agcu_phase_diagram_validated.png){ width=80% }

# Acknowledgements

The authors acknowledge the openCALPHAD project for providing the reference implementation and TDB files used for validation.

# References
