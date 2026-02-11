# DifferentiationInterface Proof of Concept
#
# This script tests whether DifferentiationInterface can replace ForwardDiff
# in OpenCALPHAD.jl with equivalent results.
#
# Usage:
#   cd OpenCALPHAD.jl
#   julia --project -e 'using Pkg; Pkg.add("DifferentiationInterface")'
#   julia --project test/di_poc.jl

using LinearAlgebra

# === Load packages ===
println("Loading packages...")
using ForwardDiff
using DifferentiationInterface
import DifferentiationInterface as DI

println("  ForwardDiff: v", pkgversion(ForwardDiff))
println("  DifferentiationInterface: v", pkgversion(DifferentiationInterface))

# === Test 1: Simple scalar derivative ===
println("\n" * "="^60)
println("Test 1: Scalar derivative (f(x) = x^3 + 2x)")
println("="^60)

f1(x) = x^3 + 2x
x_test = 3.0

# ForwardDiff direct
df_fd = ForwardDiff.derivative(f1, x_test)

# DifferentiationInterface with ForwardDiff backend
df_di = DI.derivative(f1, DI.AutoForwardDiff(), x_test)

# Analytical: f'(x) = 3x^2 + 2 => f'(3) = 29
df_analytical = 3 * x_test^2 + 2

println("  x = $x_test")
println("  ForwardDiff.derivative:  $df_fd")
println("  DI.derivative:           $df_di")
println("  Analytical:              $df_analytical")
println("  Match: $(isapprox(df_fd, df_di, rtol=1e-10) ? "✓" : "✗")")

# === Test 2: Gradient of vector function ===
println("\n" * "="^60)
println("Test 2: Gradient (f(x) = sum(x.^2) + prod(x))")
println("="^60)

f2(x) = sum(x .^ 2) + prod(x)
x_vec = [1.0, 2.0, 3.0]

# ForwardDiff direct
grad_fd = ForwardDiff.gradient(f2, x_vec)

# DifferentiationInterface
grad_di = DI.gradient(f2, DI.AutoForwardDiff(), x_vec)

println("  x = $x_vec")
println("  ForwardDiff.gradient: $grad_fd")
println("  DI.gradient:          $grad_di")
println("  Match: $(isapprox(grad_fd, grad_di, rtol=1e-10) ? "✓" : "✗")")

# === Test 3: Hessian ===
println("\n" * "="^60)
println("Test 3: Hessian (same function)")
println("="^60)

# ForwardDiff direct
hess_fd = ForwardDiff.hessian(f2, x_vec)

# DifferentiationInterface
hess_di = DI.hessian(f2, DI.AutoForwardDiff(), x_vec)

println("  ForwardDiff.hessian:")
display(hess_fd)
println("\n  DI.hessian:")
display(hess_di)
println("\n  Match: $(isapprox(hess_fd, hess_di, rtol=1e-10) ? "✓" : "✗")")

# === Test 4: Second derivative (derivative of derivative) ===
println("\n" * "="^60)
println("Test 4: Second derivative (f''(x) for f(x) = x^3 + 2x)")
println("="^60)

# ForwardDiff: derivative of derivative
df_fd_func(x) = ForwardDiff.derivative(f1, x)
d2f_fd = ForwardDiff.derivative(df_fd_func, x_test)

# DifferentiationInterface: second_derivative
d2f_di = DI.second_derivative(f1, DI.AutoForwardDiff(), x_test)

# Analytical: f''(x) = 6x => f''(3) = 18
d2f_analytical = 6 * x_test

println("  x = $x_test")
println("  ForwardDiff (nested):    $d2f_fd")
println("  DI.second_derivative:    $d2f_di")
println("  Analytical:              $d2f_analytical")
println("  Match: $(isapprox(d2f_fd, d2f_di, rtol=1e-10) ? "✓" : "✗")")

# === Test 5: Anonymous function with closure (like in differentiation.jl) ===
println("\n" * "="^60)
println("Test 5: Anonymous function with closure")
println("="^60)

# Simulate the pattern in differentiation.jl
T = 1000.0
P = 1e5
y_size = (1, 2)
y_flat = [0.3, 0.7]

# Mock Gibbs energy calculation
function mock_gibbs_energy(T, y_mat, P)
    # G = RT * sum(y * log(y)) + L * y[1] * y[2]
    R = 8.314
    L = 20000.0
    G_ideal = R * T * sum(yi > 0 ? yi * log(yi) : 0.0 for yi in y_mat)
    G_excess = L * y_mat[1] * y_mat[2]
    return G_ideal + G_excess
end

# ForwardDiff pattern (current code)
grad_fd_closure = ForwardDiff.gradient(y_flat) do yv
    y_mat = reshape(yv, y_size)
    mock_gibbs_energy(T, y_mat, P)
end

# DifferentiationInterface pattern (proposed)
G_func(yv) = mock_gibbs_energy(T, reshape(yv, y_size), P)
grad_di_closure = DI.gradient(G_func, DI.AutoForwardDiff(), y_flat)

println("  y_flat = $y_flat")
println("  ForwardDiff.gradient: $grad_fd_closure")
println("  DI.gradient:          $grad_di_closure")
println("  Match: $(isapprox(grad_fd_closure, grad_di_closure, rtol=1e-10) ? "✓" : "✗")")

# === Test 6: Hessian with closure ===
println("\n" * "="^60)
println("Test 6: Hessian with closure")
println("="^60)

hess_fd_closure = ForwardDiff.hessian(y_flat) do yv
    y_mat = reshape(yv, y_size)
    mock_gibbs_energy(T, y_mat, P)
end

hess_di_closure = DI.hessian(G_func, DI.AutoForwardDiff(), y_flat)

println("  ForwardDiff.hessian:")
display(hess_fd_closure)
println("\n  DI.hessian:")
display(hess_di_closure)
println("\n  Match: $(isapprox(hess_fd_closure, hess_di_closure, rtol=1e-10) ? "✓" : "✗")")

# === Summary ===
println("\n" * "="^60)
println("Summary")
println("="^60)

all_pass = all([
    isapprox(df_fd, df_di, rtol = 1e-10),
    isapprox(grad_fd, grad_di, rtol = 1e-10),
    isapprox(hess_fd, hess_di, rtol = 1e-10),
    isapprox(d2f_fd, d2f_di, rtol = 1e-10),
    isapprox(grad_fd_closure, grad_di_closure, rtol = 1e-10),
    isapprox(hess_fd_closure, hess_di_closure, rtol = 1e-10),
])

if all_pass
    println("✓ All tests passed!")
    println("  DifferentiationInterface can replace ForwardDiff.")
else
    println("✗ Some tests failed!")
end

# =============================================================================
# Part 2: Tests with actual OpenCALPHAD.jl
# =============================================================================

println("\n" * "="^60)
println("Part 2: OpenCALPHAD.jl Integration Tests")
println("="^60)

using OpenCALPHAD

# --- Helper: Create test phase and database (from runtests.jl) ---
function create_test_phase_and_db()
    db = Database()

    # Add elements
    push!(db.elements, Element("CR", "BCC_A2", 51.996, 4050.0, 23.543))
    push!(db.elements, Element("FE", "BCC_A2", 55.845, 4489.0, 27.28))

    # Create BCC phase with single sublattice (CR,FE)
    bcc = Phase("BCC_A2", 'S', [1.0])
    push!(bcc.constituents, ["CR", "FE"])
    push!(db.phases, bcc)

    # Add Gibbs energy parameters
    # G(BCC,CR) = GHSERCR
    push!(
        bcc.parameters,
        Parameter(
            "BCC_A2", 'G', [["CR"]], 0,
            [298.15, 6000.0],
            ["-8856.94+157.48*T-26.908*T*LN(T)"],
        ),
    )
    # G(BCC,FE) = GHSERFE
    push!(
        bcc.parameters,
        Parameter(
            "BCC_A2", 'G', [["FE"]], 0,
            [298.15, 6000.0],
            ["1225.7+124.134*T-23.5143*T*LN(T)"],
        ),
    )
    # L(BCC,CR,FE;0) - interaction parameter
    push!(
        bcc.parameters,
        Parameter(
            "BCC_A2", 'L', [["CR", "FE"]], 0,
            [298.15, 6000.0],
            ["20500-9.68*T"],
        ),
    )

    return bcc, db
end

# === Test 7: Gradient with actual calculate_gibbs_energy ===
println("\n" * "="^60)
println("Test 7: Gradient with calculate_gibbs_energy (Cr-Fe BCC)")
println("="^60)

bcc, db = create_test_phase_and_db()
T_test = 1000.0
P_test = 1e5
y_test = [0.3 0.7]  # 30% CR, 70% FE
y_flat_test = vec(y_test)
y_size_test = size(y_test)

# ForwardDiff (current implementation pattern)
grad_fd_real = ForwardDiff.gradient(y_flat_test) do yv
    y_mat = reshape(yv, y_size_test)
    calculate_gibbs_energy(bcc, T_test, y_mat, db; P = P_test)
end

# DifferentiationInterface (proposed pattern)
G_real_func(yv) =
    calculate_gibbs_energy(bcc, T_test, reshape(yv, y_size_test), db; P = P_test)
grad_di_real = DI.gradient(G_real_func, DI.AutoForwardDiff(), y_flat_test)

println("  Phase: BCC_A2 (Cr-Fe)")
println("  T = $T_test K, y = $y_flat_test")
println("  ForwardDiff.gradient: $grad_fd_real")
println("  DI.gradient:          $grad_di_real")
println("  Match: $(isapprox(grad_fd_real, grad_di_real, rtol=1e-10) ? "✓" : "✗")")

test7_pass = isapprox(grad_fd_real, grad_di_real, rtol = 1e-10)

# === Test 8: Hessian with actual calculate_gibbs_energy ===
println("\n" * "="^60)
println("Test 8: Hessian with calculate_gibbs_energy (Cr-Fe BCC)")
println("="^60)

# ForwardDiff
hess_fd_real = ForwardDiff.hessian(y_flat_test) do yv
    y_mat = reshape(yv, y_size_test)
    calculate_gibbs_energy(bcc, T_test, y_mat, db; P = P_test)
end

# DifferentiationInterface
hess_di_real = DI.hessian(G_real_func, DI.AutoForwardDiff(), y_flat_test)

println("  ForwardDiff.hessian:")
display(hess_fd_real)
println("\n  DI.hessian:")
display(hess_di_real)
println("\n  Match: $(isapprox(hess_fd_real, hess_di_real, rtol=1e-10) ? "✓" : "✗")")

test8_pass = isapprox(hess_fd_real, hess_di_real, rtol = 1e-10)

# === Test 9: Scalar derivative (like two_phase.jl pattern) ===
println("\n" * "="^60)
println("Test 9: Scalar derivative (two_phase.jl pattern)")
println("="^60)

# G as function of x (first component mole fraction)
function G_at_x(x)
    n_const = 2
    n_sub = 1
    y = zeros(typeof(x), n_sub, n_const)
    y[1, 1] = x
    y[1, 2] = one(x) - x
    return calculate_gibbs_energy(bcc, T_test, y, db; P = P_test)
end

x_scalar = 0.3

# ForwardDiff (current pattern in two_phase.jl)
dG_fd = ForwardDiff.derivative(G_at_x, x_scalar)

# DifferentiationInterface
dG_di = DI.derivative(G_at_x, DI.AutoForwardDiff(), x_scalar)

println("  x = $x_scalar (mole fraction of CR)")
println("  ForwardDiff.derivative: $dG_fd")
println("  DI.derivative:          $dG_di")
println("  Match: $(isapprox(dG_fd, dG_di, rtol=1e-10) ? "✓" : "✗")")

test9_pass = isapprox(dG_fd, dG_di, rtol = 1e-10)

# === Test 10: Second derivative (d²G/dx²) ===
println("\n" * "="^60)
println("Test 10: Second derivative d²G/dx² (diffusion potential)")
println("="^60)

# ForwardDiff (current pattern)
dG_func(x) = ForwardDiff.derivative(G_at_x, x)
d2G_fd = ForwardDiff.derivative(dG_func, x_scalar)

# DifferentiationInterface
d2G_di = DI.second_derivative(G_at_x, DI.AutoForwardDiff(), x_scalar)

println("  x = $x_scalar")
println("  ForwardDiff (nested):    $d2G_fd")
println("  DI.second_derivative:    $d2G_di")
println("  Match: $(isapprox(d2G_fd, d2G_di, rtol=1e-10) ? "✓" : "✗")")

test10_pass = isapprox(d2G_fd, d2G_di, rtol = 1e-10)

# === Test 11: Compare with existing compute_gradient ===
println("\n" * "="^60)
println("Test 11: Compare with existing compute_gradient/compute_hessian")
println("="^60)

# Use existing OpenCALPHAD functions
grad_existing = compute_gradient(bcc, T_test, y_test, db, ForwardDiffMethod(); P = P_test)
hess_existing = compute_hessian(bcc, T_test, y_test, db, ForwardDiffMethod(); P = P_test)

println("  compute_gradient (existing): $grad_existing")
println("  DI.gradient:                 $grad_di_real")
println(
    "  Gradient match: $(isapprox(grad_existing, grad_di_real, rtol=1e-10) ? "✓" : "✗")",
)

println("\n  compute_hessian (existing):")
display(hess_existing)
println("\n  DI.hessian:")
display(hess_di_real)
println(
    "\n  Hessian match: $(isapprox(hess_existing, hess_di_real, rtol=1e-10) ? "✓" : "✗")",
)

test11_grad_pass = isapprox(grad_existing, grad_di_real, rtol = 1e-10)
test11_hess_pass = isapprox(hess_existing, hess_di_real, rtol = 1e-10)

# === Final Summary ===
println("\n" * "="^60)
println("Final Summary")
println("="^60)

all_basic_pass = all([
    isapprox(df_fd, df_di, rtol = 1e-10),
    isapprox(grad_fd, grad_di, rtol = 1e-10),
    isapprox(hess_fd, hess_di, rtol = 1e-10),
    isapprox(d2f_fd, d2f_di, rtol = 1e-10),
    isapprox(grad_fd_closure, grad_di_closure, rtol = 1e-10),
    isapprox(hess_fd_closure, hess_di_closure, rtol = 1e-10),
])

all_opencalphad_pass = all([
    test7_pass,
    test8_pass,
    test9_pass,
    test10_pass,
    test11_grad_pass,
    test11_hess_pass,
])

println("Part 1 (Basic tests):       $(all_basic_pass ? "✓ All passed" : "✗ Some failed")")
println(
    "Part 2 (OpenCALPHAD tests): $(all_opencalphad_pass ? "✓ All passed" : "✗ Some failed")",
)

if all_basic_pass && all_opencalphad_pass
    println("\n✓ All tests passed!")
    println("  DifferentiationInterface can safely replace ForwardDiff in OpenCALPHAD.jl")
else
    println("\n✗ Some tests failed!")
end

println("\nAPI comparison:")
println("  ForwardDiff.derivative(f, x)  →  DI.derivative(f, DI.AutoForwardDiff(), x)")
println("  ForwardDiff.gradient(f, x)    →  DI.gradient(f, DI.AutoForwardDiff(), x)")
println("  ForwardDiff.hessian(f, x)     →  DI.hessian(f, DI.AutoForwardDiff(), x)")

# =============================================================================
# Part 3: DIMethod vs ForwardDiffMethod (after DIMethod implementation)
# =============================================================================

println("\n" * "="^60)
println("Part 3: DIMethod vs ForwardDiffMethod")
println("="^60)

# Check if DIMethod is defined
di_method_available = isdefined(OpenCALPHAD, :DIMethod)

if !di_method_available
    println("  DIMethod not yet implemented - skipping Tests 12-13")
    println("  (Run this script again after implementing DIMethod)")
    test12_pass = missing
    test13_pass = missing
else
    # === Test 12: DIMethod gradient vs ForwardDiffMethod gradient ===
    println("\n" * "="^60)
    println("Test 12: compute_gradient - DIMethod vs ForwardDiffMethod")
    println("="^60)

    grad_fdm = compute_gradient(bcc, T_test, y_test, db, ForwardDiffMethod(); P = P_test)
    grad_dim = compute_gradient(bcc, T_test, y_test, db, DIMethod(); P = P_test)

    println("  Phase: BCC_A2 (Cr-Fe), T = $T_test K")
    println("  ForwardDiffMethod: $grad_fdm")
    println("  DIMethod:          $grad_dim")
    println("  Match: $(isapprox(grad_fdm, grad_dim, rtol=1e-10) ? "✓" : "✗")")

    test12_pass = isapprox(grad_fdm, grad_dim, rtol = 1e-10)

    # === Test 13: DIMethod hessian vs ForwardDiffMethod hessian ===
    println("\n" * "="^60)
    println("Test 13: compute_hessian - DIMethod vs ForwardDiffMethod")
    println("="^60)

    hess_fdm = compute_hessian(bcc, T_test, y_test, db, ForwardDiffMethod(); P = P_test)
    hess_dim = compute_hessian(bcc, T_test, y_test, db, DIMethod(); P = P_test)

    println("  ForwardDiffMethod:")
    display(hess_fdm)
    println("\n  DIMethod:")
    display(hess_dim)
    println("\n  Match: $(isapprox(hess_fdm, hess_dim, rtol=1e-10) ? "✓" : "✗")")

    test13_pass = isapprox(hess_fdm, hess_dim, rtol = 1e-10)
end

# =============================================================================
# Part 4: phase_field.jl functions (Test 14-15)
# =============================================================================

println("\n" * "="^60)
println("Part 4: phase_field.jl Functions")
println("="^60)

# === Test 14: diffusion_potential comparison ===
println("\n" * "="^60)
println("Test 14: diffusion_potential (d²G/dx²)")
println("="^60)

# Current implementation uses ForwardDiff.derivative nested
# We compare with DI.second_derivative

x_pf = 0.3

# Simulate current diffusion_potential implementation
function diffusion_potential_fd(phase, T, x, db; P = 1e5)
    function G_func(x_val)
        y = zeros(typeof(x_val), 1, 2)
        y[1, 1] = x_val
        y[1, 2] = one(x_val) - x_val
        return calculate_gibbs_energy(phase, T, y, db; P = P)
    end
    dGdx(x_val) = ForwardDiff.derivative(G_func, x_val)
    return ForwardDiff.derivative(dGdx, x)
end

# Proposed DI implementation
function diffusion_potential_di(phase, T, x, db; P = 1e5)
    function G_func(x_val)
        y = zeros(typeof(x_val), 1, 2)
        y[1, 1] = x_val
        y[1, 2] = one(x_val) - x_val
        return calculate_gibbs_energy(phase, T, y, db; P = P)
    end
    return DI.second_derivative(G_func, DI.AutoForwardDiff(), x)
end

d2G_fd_pf = diffusion_potential_fd(bcc, T_test, x_pf, db; P = P_test)
d2G_di_pf = diffusion_potential_di(bcc, T_test, x_pf, db; P = P_test)

# Also compare with existing diffusion_potential function
d2G_existing = diffusion_potential(bcc, T_test, x_pf, db; P = P_test)

println("  Phase: BCC_A2 (Cr-Fe), T = $T_test K, x = $x_pf")
println("  ForwardDiff (nested):      $d2G_fd_pf")
println("  DI.second_derivative:      $d2G_di_pf")
println("  diffusion_potential():     $d2G_existing")
println("  FD vs DI match: $(isapprox(d2G_fd_pf, d2G_di_pf, rtol=1e-10) ? "✓" : "✗")")
println("  Existing match: $(isapprox(d2G_existing, d2G_di_pf, rtol=1e-10) ? "✓" : "✗")")

test14_pass =
    isapprox(d2G_fd_pf, d2G_di_pf, rtol = 1e-10) &&
    isapprox(d2G_existing, d2G_di_pf, rtol = 1e-10)

# === Test 15: chemical_potential derivative pattern ===
println("\n" * "="^60)
println("Test 15: chemical_potential derivative pattern (dG/dx)")
println("="^60)

# Simulate the derivative pattern used in chemical_potential
function chem_pot_derivative_fd(phase, T, x, db; P = 1e5)
    function G_func(x_val)
        y = zeros(typeof(x_val), 1, 2)
        y[1, 1] = x_val
        y[1, 2] = one(x_val) - x_val
        return calculate_gibbs_energy(phase, T, y, db; P = P)
    end
    G = G_func(x)
    dGdx = ForwardDiff.derivative(G_func, x)
    μ₁ = G - x * dGdx
    μ₂ = G + (1 - x) * dGdx
    return (μ₁, μ₂)
end

function chem_pot_derivative_di(phase, T, x, db; P = 1e5)
    function G_func(x_val)
        y = zeros(typeof(x_val), 1, 2)
        y[1, 1] = x_val
        y[1, 2] = one(x_val) - x_val
        return calculate_gibbs_energy(phase, T, y, db; P = P)
    end
    G = G_func(x)
    dGdx = DI.derivative(G_func, DI.AutoForwardDiff(), x)
    μ₁ = G - x * dGdx
    μ₂ = G + (1 - x) * dGdx
    return (μ₁, μ₂)
end

mu_fd = chem_pot_derivative_fd(bcc, T_test, x_pf, db; P = P_test)
mu_di = chem_pot_derivative_di(bcc, T_test, x_pf, db; P = P_test)

# Compare with existing chemical_potential (db version)
mu_existing = chemical_potential(bcc, T_test, x_pf, db; P = P_test)

println("  Phase: BCC_A2 (Cr-Fe), T = $T_test K, x = $x_pf")
println("  ForwardDiff: μ₁=$(mu_fd[1]), μ₂=$(mu_fd[2])")
println("  DI:          μ₁=$(mu_di[1]), μ₂=$(mu_di[2])")
println("  Existing:    μ₁=$(mu_existing[1]), μ₂=$(mu_existing[2])")
mu_fd_di_match =
    isapprox(mu_fd[1], mu_di[1], rtol = 1e-10) && isapprox(mu_fd[2], mu_di[2], rtol = 1e-10)
mu_existing_di_match =
    isapprox(mu_existing[1], mu_di[1], rtol = 1e-10) &&
    isapprox(mu_existing[2], mu_di[2], rtol = 1e-10)

println("  FD vs DI match: $(mu_fd_di_match ? "✓" : "✗")")
println("  Existing match: $(mu_existing_di_match ? "✓" : "✗ (different implementation)")")

# Note: Existing chemical_potential uses gibbs_energy_with_derivatives which may differ
# The important test is FD vs DI match
test15_pass = mu_fd_di_match

all_phasefield_pass = test14_pass && test15_pass

# === Grand Summary ===
println("\n" * "="^60)
println("Grand Summary")
println("="^60)

println("Part 1 (Basic DI tests):    $(all_basic_pass ? "✓ All passed" : "✗ Some failed")")
println(
    "Part 2 (OpenCALPHAD + DI):  $(all_opencalphad_pass ? "✓ All passed" : "✗ Some failed")",
)

if di_method_available
    all_dimethod_pass = test12_pass && test13_pass
    println(
        "Part 3 (DIMethod dispatch): $(all_dimethod_pass ? "✓ All passed" : "✗ Some failed")",
    )
else
    println("Part 3 (DIMethod dispatch): ⏳ Pending (DIMethod not implemented)")
end

println(
    "Part 4 (phase_field.jl):    $(all_phasefield_pass ? "✓ All passed" : "✗ Some failed")",
)

if di_method_available
    all_tests_pass =
        all_basic_pass && all_opencalphad_pass && all_dimethod_pass && all_phasefield_pass
else
    all_tests_pass = all_basic_pass && all_opencalphad_pass && all_phasefield_pass
end

if all_tests_pass
    if di_method_available
        println("\n✓ All tests passed! DIMethod implementation verified.")
    else
        println("\n✓ All DI API tests passed!")
        println("Next step: Implement DIMethod in src/Types/functions.jl")
    end
else
    println("\n✗ Some tests failed!")
end
