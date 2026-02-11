# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Magnetic parameter optimization
#
# Magnetic Parameter Optimization Example
#
# Demonstrates that ForwardDiff can differentiate through
# conditional branching (Inden-Hillert model) for parameter optimization.
#
# We optimize Tc (Curie temperature) and beta (magnetic moment)
# to fit target magnetic contribution data.
#
# Usage:
#   julia --project examples/624_magnetic_optimization.jl

using OpenCALPHAD
using Optim
using ForwardDiff
using Plots
using Printf

println("=" ^ 70)
println("Magnetic Parameter Optimization (Tc, β)")
println("=" ^ 70)

const R = 8.314462618

# =============================================================================
# Inden-Hillert-Jarl Model (same as before)
# =============================================================================

function inden_hillert_g(tau, p)
    A = 518/1125 + (11692/15975) * (1/p - 1)
    if tau < 1
        term1 = 79 / (140 * p)
        term2 = (474/497) * (1/p - 1) * (tau^3/6 + tau^9/135 + tau^15/600)
        g = 1 - (term1 * tau^(-1) + term2) / A
    else
        g = -(tau^(-5)/10 + tau^(-15)/315 + tau^(-25)/1500) / A
    end
    return g
end

function G_magnetic(T, Tc, beta, p)
    tau = T / Tc
    g = inden_hillert_g(tau, p)
    return R * T * log(beta + 1) * g
end

# =============================================================================
# Target data (generated from "true" parameters)
# =============================================================================

# True parameters (we pretend we don't know these)
const Tc_true = 1043.0
const beta_true = 2.22
const p_bcc = 0.40

# Generate target data at various temperatures
const T_data = [400.0, 600.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1400.0]
const G_target = [G_magnetic(T, Tc_true, beta_true, p_bcc) for T in T_data]

println("\nTarget data (from Tc=$Tc_true K, β=$beta_true):")
for (T, G) in zip(T_data, G_target)
    @printf("  T = %6.0f K: G_mag = %8.2f J/mol\n", T, G)
end

# =============================================================================
# Loss function
# =============================================================================

function loss(params)
    Tc, beta = params

    # Bounds check (Tc and beta must be positive)
    if Tc <= 0 || beta <= 0
        return 1e10
    end

    total_error = 0.0
    for (T, G_tgt) in zip(T_data, G_target)
        G_calc = G_magnetic(T, Tc, beta, p_bcc)
        total_error += (G_calc - G_tgt)^2
    end
    return total_error
end

# =============================================================================
# Test gradient computation through conditional
# =============================================================================

println("\n" * "-" ^ 70)
println("Testing gradient computation through conditional branching")
println("-" ^ 70)

# Initial guess (perturbed from true values)
Tc_init = 1100.0   # +5.5% error
beta_init = 2.5    # +12.6% error
initial_params = [Tc_init, beta_init]

println("\nInitial parameters:")
@printf("  Tc = %.1f K (true: %.1f K, error: %.1f%%)\n",
    Tc_init, Tc_true, (Tc_init - Tc_true)/Tc_true * 100)
@printf("  β = %.2f (true: %.2f, error: %.1f%%)\n",
    beta_init, beta_true, (beta_init - beta_true)/beta_true * 100)

# Compute gradient using ForwardDiff
grad = ForwardDiff.gradient(loss, initial_params)
println("\nGradient at initial point (ForwardDiff):")
@printf("  ∂loss/∂Tc = %.2e\n", grad[1])
@printf("  ∂loss/∂β  = %.2e\n", grad[2])

# Verify with finite differences
h = 1e-5
grad_fd = [
    (loss([Tc_init + h, beta_init]) - loss([Tc_init - h, beta_init])) / (2h),
    (loss([Tc_init, beta_init + h]) - loss([Tc_init, beta_init - h])) / (2h),
]
println("\nGradient at initial point (Finite Difference):")
@printf("  ∂loss/∂Tc = %.2e\n", grad_fd[1])
@printf("  ∂loss/∂β  = %.2e\n", grad_fd[2])

println("\nRelative difference:")
@printf("  Tc: %.2e\n", abs(grad[1] - grad_fd[1]) / abs(grad_fd[1]))
@printf("  β:  %.2e\n", abs(grad[2] - grad_fd[2]) / abs(grad_fd[2]))

# =============================================================================
# Optimization with BFGS (gradient-based)
# =============================================================================

println("\n" * "-" ^ 70)
println("Optimization with BFGS (gradient-based, uses ForwardDiff)")
println("-" ^ 70)

# Define gradient function using ForwardDiff
function loss_grad!(G, params)
    G .= ForwardDiff.gradient(loss, params)
end

result_bfgs = optimize(loss, loss_grad!, initial_params, BFGS(),
    Optim.Options(show_trace = false, iterations = 100))

Tc_opt, beta_opt = Optim.minimizer(result_bfgs)

println("\nBFGS Result:")
@printf("  Tc = %.2f K (true: %.1f K, error: %.4f%%)\n",
    Tc_opt, Tc_true, abs(Tc_opt - Tc_true)/Tc_true * 100)
@printf("  β = %.4f (true: %.2f, error: %.4f%%)\n",
    beta_opt, beta_true, abs(beta_opt - beta_true)/beta_true * 100)
@printf("  Iterations: %d\n", Optim.iterations(result_bfgs))
@printf("  Final loss: %.2e\n", Optim.minimum(result_bfgs))

# =============================================================================
# Comparison with Nelder-Mead (derivative-free)
# =============================================================================

println("\n" * "-" ^ 70)
println("Comparison with Nelder-Mead (derivative-free)")
println("-" ^ 70)

result_nm = optimize(loss, initial_params, NelderMead(),
    Optim.Options(show_trace = false, iterations = 200))

Tc_nm, beta_nm = Optim.minimizer(result_nm)

println("\nNelder-Mead Result:")
@printf("  Tc = %.2f K (error: %.4f%%)\n",
    Tc_nm, abs(Tc_nm - Tc_true)/Tc_true * 100)
@printf("  β = %.4f (error: %.4f%%)\n",
    beta_nm, abs(beta_nm - beta_true)/beta_true * 100)
@printf("  Iterations: %d\n", Optim.iterations(result_nm))
@printf("  Final loss: %.2e\n", Optim.minimum(result_nm))

# =============================================================================
# Plot results
# =============================================================================

println("\nGenerating comparison plot...")

T_plot = 300.0:10.0:1600.0
G_true = [G_magnetic(T, Tc_true, beta_true, p_bcc) for T in T_plot]
G_init = [G_magnetic(T, Tc_init, beta_init, p_bcc) for T in T_plot]
G_bfgs = [G_magnetic(T, Tc_opt, beta_opt, p_bcc) for T in T_plot]

p = plot(T_plot, G_true,
    label = "Target (Tc=1043, β=2.22)",
    xlabel = "Temperature [K]",
    ylabel = "G_magnetic [J/mol]",
    title = "Magnetic Parameter Optimization",
    linewidth = 2,
    color = :black,
    linestyle = :solid,
    size = (700, 500),
)

plot!(p, T_plot, G_init,
    label = "Initial (Tc=1100, β=2.50)",
    linewidth = 2,
    color = :red,
    linestyle = :dash,
)

plot!(p, T_plot, G_bfgs,
    label = @sprintf("BFGS result (Tc=%.1f, β=%.2f)", Tc_opt, beta_opt),
    linewidth = 2,
    color = :blue,
    linestyle = :solid,
)

scatter!(p, T_data, G_target,
    label = "Data points",
    markersize = 8,
    color = :black,
    markershape = :circle,
)

vline!(p, [Tc_true], linestyle = :dot, color = :gray, label = "")
annotate!(p, Tc_true + 30, -3000, text("Tc", 9, :left))

output_path = joinpath(@__DIR__, "624_magnetic_optimization.png")
savefig(p, output_path)
println("Saved to: $output_path")

# =============================================================================
# Summary
# =============================================================================

println("\n" * "=" ^ 70)
println("Summary")
println("=" ^ 70)

println(
    """

Key Result: ForwardDiff successfully computes gradients through
the Inden-Hillert conditional branching (τ < 1 vs τ ≥ 1).

Method Comparison:
                    BFGS (autodiff)    Nelder-Mead
  Iterations:       $(lpad(Optim.iterations(result_bfgs), 10))       $(lpad(Optim.iterations(result_nm), 10))
  Final loss:       $(lpad(@sprintf("%.2e", Optim.minimum(result_bfgs)), 10))       $(lpad(@sprintf("%.2e", Optim.minimum(result_nm)), 10))

BFGS with autodiff converges faster because it uses exact gradients
computed by ForwardDiff through the conditional branching.

This capability is not possible with symbolic differentiation,
which cannot handle if/else control flow.
""",
)

display(p)
