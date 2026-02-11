# Reference Tests for OpenCALPHAD.jl

Last-Modified: 2026-02-11

This directory contains reference tests that verify OpenCALPHAD.jl produces
the same results as the original openCALPHAD (Fortran).

## Overview

The reference tests compare OpenCALPHAD.jl calculations against pre-computed
values from openCALPHAD (Fortran). This ensures numerical equivalence and
serves as a regression test for the Julia implementation.

### What is Tested

| Category | What is Compared | Tolerance | Status |
|----------|------------------|-----------|--------|
| TDB parsing | Phase names, constituents, parameters | exact | Tested |
| Gibbs energy | G(T, x) | ±10 J/mol | Tested |
| Chemical potential | μ(T, x) | ±10 J/mol | Tested |
| Equilibrium | Phase fractions, compositions | ±0.01 | Tested |
| STEP | T-x diagram points | ±1 K, ±0.01 | Tested |

### Test Systems

| System | TDB File | Phases | Description |
|--------|----------|--------|-------------|
| Ag-Cu | agcu.TDB | FCC_A1, LIQUID | Miscibility gap equilibrium |

## Quick Start

### Running Tests (No openCALPHAD needed)

```bash
# Navigate to the OpenCALPHAD.jl directory
cd OpenCALPHAD.jl

# Lightweight tests only (CI, ~10 seconds)
julia --project=. -e 'using Pkg; Pkg.test()'

# Numerical verification tests (~20 minutes)
julia --project=. reftest/run_numerical_tests.jl

# Fortran reference comparison tests
julia --project=. reftest/run_reftest.jl
```

### Test Organization

| Script | Tests | Time | Purpose |
|--------|-------|------|---------|
| `Pkg.test()` | 97 | ~10s | Type construction, error fallbacks, solver basics |
| `run_numerical_tests.jl` | 109 | ~20m | Equilibrium, STEP, MAP, plotting, multicomponent |
| `run_reftest.jl` | — | ~5m | Fortran reference value comparison |

### Regenerating Reference Data (Requires openCALPHAD)

If you want to verify the reference values yourself:

```bash
# 1. Clone and build openCALPHAD
git clone https://github.com/sundmanbo/opencalphad
cd opencalphad
make
export PATH=$PATH:$(pwd)/bin

# 2. Navigate to reftest directory
cd /path/to/OpenCALPHAD.jl/reftest

# 3. Generate reference data
./generate_reference.sh agcu_equil_1000K

# 4. Review the log output
cat reference/agcu_equil_1000K.log
```

## Reference Implementation

The reference values are generated using openCALPHAD (Fortran).

### Package Provenance

| Item | Value |
|------|-------|
| Repository | https://github.com/sundmanbo/opencalphad |
| Commit | `e857325f5d38e32e41513186ddf1d8ae4b2e0349` |
| Date | 2026-01-01 |
| Version | OpenCALPHAD 6.x |
| Platform | macOS arm64 (Apple Silicon) |

### Local Clone Path

```
../opencalphad/
```

## Directory Structure

```
reftest/
├── README.md                 # This file
├── run_reftest.jl            # Fortran reference comparison tests
├── run_numerical_tests.jl    # Numerical verification tests (separated from runtests.jl)
├── generate_reference.sh     # Shell script to generate reference data
├── macros/                   # OpenCALPHAD macro files (.OCM)
│   └── agcu_equil_1000K.OCM
├── reference/                # Pre-computed reference data (.json, .log)
│   ├── agcu_equil_1000K.json
│   └── agcu_equil_1000K.log
└── tdb/                      # Thermodynamic database files
    ├── agcu.TDB
    └── steel1.TDB
```

## Reference Data Format

Each JSON file contains:

```json
{
  "metadata": {
    "generator": "openCALPHAD",
    "version": "6.100",
    "platform": "macOS arm64",
    "generated": "2026-01-16T20:00:00+09:00",
    "macro": "macros/agcu_equil_1000K.OCM",
    "tdb": "tdb/agcu.TDB",
    "reproduction": "cd reftest && ./generate_reference.sh agcu_equil_1000K"
  },
  "conditions": {
    "T": 1000.0,
    "P": 100000.0,
    "x_Cu": 0.3
  },
  "results": {
    "phases": [...],
    "equilibrium": {...}
  }
}
```

## Available Test Cases

| Test Name | System | Conditions | Description |
|-----------|--------|------------|-------------|
| `agcu_equil_1000K` | Ag-Cu | T=1000K, x(Cu)=0.3 | Miscibility gap equilibrium |

## Tolerance Rationale

The tolerances are chosen based on typical CALPHAD precision:

| Property | Tolerance | Rationale |
|----------|-----------|-----------|
| Gibbs energy | ±10 J/mol | Typical assessment uncertainty |
| Phase amounts | ±0.01 mol | 1% precision |
| Compositions | ±0.01 | Mole fraction precision |

These tolerances are more relaxed than numerical precision (which would be
~1e-10) because CALPHAD calculations involve:
1. Iterative solvers with different convergence criteria
2. Slightly different numerical methods between implementations
3. Round-off accumulation in complex expressions

## Adding New Test Cases

### Step 1: Create OCM Macro

Create a new macro file in `macros/`:

```ocm
$ macros/new_test.OCM
$ Description of test
new Y
read tdb ./tdb/system.TDB
set cond t=1000 p=1e5 x(B)=0.3
c e
list eq 1
```

### Step 2: Add TDB File

If needed, add the TDB file to `tdb/`.

### Step 3: Generate Reference Data

```bash
./generate_reference.sh new_test
```

### Step 4: Parse and Create JSON

Review `reference/new_test.log` and create `reference/new_test.json` with
the structured data.

### Step 5: Add Test to run_reftest.jl

Add a new test case in the appropriate @testset block.

## Comparison with Other Implementations

### vs pycalphad

OpenCALPHAD.jl aims to match openCALPHAD (Fortran) rather than pycalphad
because:
1. Same TDB parsing logic as the original
2. Same numerical methods (CEF model)
3. Direct comparison with reference implementation

### vs Thermo-Calc

Commercial software with different internal algorithms. Not used for
reference testing.

## Troubleshooting

### Test Fails with Large Error

1. Check if the TDB file matches the reference
2. Verify temperature and composition conditions
3. Check if the phase model is fully implemented

### Cannot Regenerate Reference Data

1. Ensure openCALPHAD is in PATH: `which oc6`
2. Check macro file syntax
3. Review error messages in the log file

## License

- **OpenCALPHAD.jl**: GPL-3.0-or-later
- **TDB files**: Various SGTE/CALPHAD databases (see file headers)
- **Reference data**: Generated outputs, no copyright

## References

1. Sundman, B., et al. "OpenCALPHAD - a free thermodynamic software."
   Integrating Materials and Manufacturing Innovation 4:1 (2015).

2. Lukas, H. L., Fries, S. G., & Sundman, B. "Computational Thermodynamics:
   The Calphad Method." Cambridge University Press, 2007.

## See Also

- [examples/901_validation_comparison.jl](../examples/901_validation_comparison.jl) - Visual comparison
- [examples/911_phase_diagram_validated.jl](../examples/911_phase_diagram_validated.jl) - Validated phase diagram
