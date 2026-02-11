# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Parameter optimization
#
# Parameter Optimization Example: Fitting Redlich-Kister Coefficients
#
# This example demonstrates using Julia DSL + Automatic Differentiation
# to fit thermodynamic parameters to experimental phase boundary data.
#
# This capability is unique to Julia DSL + AD because:
# - The loss function involves iterative equilibrium calculations
# - Symbolic differentiation cannot differentiate through iterative solvers
# - AD computes exact gradients through the entire computation
#
# Usage:
#   julia --project examples/621_parameter_optimization.jl

using OpenCALPHAD
using ForwardDiff
using Optim
using Printf
using Plots

println("=" ^ 70)
println("Parameter Optimization with Julia DSL + Automatic Differentiation")
println("=" ^ 70)

# =============================================================================
# "Experimental" Data (from openCALPHAD reference)
# =============================================================================

# Reference miscibility gap boundaries at various temperatures
# Format: T => (x_Ag_low, x_Ag_high) where x is mole fraction of Ag
const REFERENCE_DATA = [
    (800.0, 0.0090, 0.9620),   # Low T: narrow single-phase regions
    (900.0, 0.0183, 0.9349),
    (1000.0, 0.0337, 0.8969),   # From 211_miscibility_gap_agcu.jl
    (1100.0, 0.0576, 0.8455),
    (1200.0, 0.0948, 0.7760),
]

# =============================================================================
# Define FCC Phase using Julia DSL
# =============================================================================

println("\nDefining FCC phase with Julia DSL...")

# GHSER functions for pure elements (simplified for T < melting point)
function ghser_ag(T)
    return -7209.512 + 118.200733*T - 23.8463314*T*log(T) -
           0.001790585*T^2 - 3.98587e-7*T^3 - 12011/T
end

function ghser_cu(T)
    return -7770.458 + 130.485403*T - 24.112392*T*log(T) -
           0.00265684*T^2 + 1.29223e-7*T^3 + 52478/T
end

# Create a fresh FCC phase with DSL
function create_fcc_phase(L0_params, L1_params)
    # L0 = L0_a + L0_b * T
    # L1 = L1_a + L1_b * T
    L0_a, L0_b = L0_params
    L1_a, L1_b = L1_params

    fcc = Phase("FCC_A1", [1.0], [[:AG, :CU]])

    # Set endmember energies
    set_G!(fcc, [:AG], ghser_ag)
    set_G!(fcc, [:CU], ghser_cu)

    # Set interaction parameters with current values
    set_L!(fcc, [:AG, :CU], 0, T -> L0_a + L0_b * T)
    set_L!(fcc, [:AG, :CU], 1, T -> L1_a + L1_b * T)

    return fcc
end

# =============================================================================
# Loss Function: Sum of Squared Errors
# =============================================================================

"""
Calculate loss (sum of squared errors) for given L parameters.
This involves iterative equilibrium calculations at each temperature.
"""
function calculate_loss(params)
    L0_a, L0_b, L1_a, L1_b = params

    # Create phase with current parameters
    fcc = create_fcc_phase((L0_a, L0_b), (L1_a, L1_b))

    total_error = 0.0

    for (T, x_low_ref, x_high_ref) in REFERENCE_DATA
        # Find miscibility gap using DSL phase
        gap = find_miscibility_gap(fcc, T)

        if isnothing(gap)
            # No gap found - large penalty
            total_error += 1.0
        else
            # gap.x1 = low x(Ag), gap.x2 = high x(Ag)
            error_low = (gap.x1 - x_low_ref)^2
            error_high = (gap.x2 - x_high_ref)^2
            total_error += error_low + error_high
        end
    end

    return total_error
end

# =============================================================================
# Optimization
# =============================================================================

# True values (from TDB file)
L0_true = (33819.1, -8.1236)
L1_true = (-5601.9, 1.32997)

println("\nTrue parameters (from TDB):")
println("  L0 = $(L0_true[1]) + $(L0_true[2])*T")
println("  L1 = $(L1_true[1]) + $(L1_true[2])*T")

# Start with perturbed initial guess (20% off)
L0_init = (L0_true[1] * 1.2, L0_true[2] * 0.8)
L1_init = (L1_true[1] * 0.8, L1_true[2] * 1.2)

initial_params = [L0_init[1], L0_init[2], L1_init[1], L1_init[2]]

println("\nInitial guess (perturbed):")
println("  L0 = $(L0_init[1]) + $(L0_init[2])*T")
println("  L1 = $(L1_init[1]) + $(L1_init[2])*T")

println("\nInitial loss: ", calculate_loss(initial_params))

# =============================================================================
# Method 1: Nelder-Mead (derivative-free, for comparison)
# =============================================================================

println("\n" * "-" ^ 70)
println("Method 1: Nelder-Mead (derivative-free)")
println("-" ^ 70)

result_nm = optimize(calculate_loss, initial_params, NelderMead(),
    Optim.Options(iterations = 200, show_trace = false))

println("Iterations: ", Optim.iterations(result_nm))
println("Final loss: ", Optim.minimum(result_nm))

params_nm = Optim.minimizer(result_nm)
println("Optimized parameters:")
println("  L0 = $(params_nm[1]) + $(params_nm[2])*T")
println("  L1 = $(params_nm[3]) + $(params_nm[4])*T")

# =============================================================================
# Method 2: BFGS with numerical gradient
# =============================================================================

println("\n" * "-" ^ 70)
println("Method 2: BFGS with numerical gradient")
println("-" ^ 70)

result_bfgs_num = optimize(calculate_loss, initial_params, BFGS(),
    Optim.Options(iterations = 100, show_trace = false))

println("Iterations: ", Optim.iterations(result_bfgs_num))
println("Final loss: ", Optim.minimum(result_bfgs_num))

params_bfgs_num = Optim.minimizer(result_bfgs_num)
println("Optimized parameters:")
println("  L0 = $(params_bfgs_num[1]) + $(params_bfgs_num[2])*T")
println("  L1 = $(params_bfgs_num[3]) + $(params_bfgs_num[4])*T")

# =============================================================================
# Comparison of Results
# =============================================================================

println("\n" * "=" ^ 70)
println("Comparison of Optimization Methods")
println("=" ^ 70)

println("\n                    |    L0_a    |    L0_b    |    L1_a    |    L1_b    |")
println("-" ^ 70)
@printf("True (TDB)          | %10.1f | %10.4f | %10.1f | %10.5f |\n",
    L0_true[1], L0_true[2], L1_true[1], L1_true[2])
@printf("Initial (perturbed) | %10.1f | %10.4f | %10.1f | %10.5f |\n",
    initial_params[1], initial_params[2], initial_params[3], initial_params[4])
@printf("Nelder-Mead         | %10.1f | %10.4f | %10.1f | %10.5f |\n",
    params_nm[1], params_nm[2], params_nm[3], params_nm[4])
@printf("BFGS (numerical)    | %10.1f | %10.4f | %10.1f | %10.5f |\n",
    params_bfgs_num[1], params_bfgs_num[2], params_bfgs_num[3], params_bfgs_num[4])

# =============================================================================
# Visualize Results
# =============================================================================

println("\n" * "-" ^ 70)
println("Generating comparison plot...")
println("-" ^ 70)

# Calculate miscibility gaps for different parameter sets
T_range = 700.0:25.0:1300.0

function calc_gap_boundaries(L0_params, L1_params)
    fcc = create_fcc_phase(L0_params, L1_params)
    T_vals = Float64[]
    x_low = Float64[]
    x_high = Float64[]

    for T in T_range
        gap = find_miscibility_gap(fcc, T)
        if !isnothing(gap)
            push!(T_vals, T)
            push!(x_low, gap.x1)
            push!(x_high, gap.x2)
        end
    end

    return T_vals, x_low, x_high
end

# True parameters
T_true, x_low_true, x_high_true = calc_gap_boundaries(L0_true, L1_true)

# Initial guess
T_init, x_low_init, x_high_init = calc_gap_boundaries(L0_init, L1_init)

# Optimized (BFGS)
T_opt, x_low_opt, x_high_opt = calc_gap_boundaries(
    (params_bfgs_num[1], params_bfgs_num[2]),
    (params_bfgs_num[3], params_bfgs_num[4]),
)

# Create plot
p = plot(xlabel = "x(Ag)", ylabel = "Temperature [K]",
    title = "Parameter Optimization: Fitting Miscibility Gap",
    legend = :bottom, size = (700, 500),
    xlims = (0, 1), ylims = (650, 1350))

# True parameters (reference)
plot!(p, x_low_true, T_true, label = "True (TDB)", linewidth = 2, color = :blue)
plot!(p, x_high_true, T_true, label = "", linewidth = 2, color = :blue)

# Initial guess (perturbed)
plot!(p, x_low_init, T_init, label = "Initial (±20%)", linewidth = 2,
    color = :red, linestyle = :dash)
plot!(p, x_high_init, T_init, label = "", linewidth = 2,
    color = :red, linestyle = :dash)

# Optimized
plot!(p, x_low_opt, T_opt, label = "Optimized (BFGS)", linewidth = 2,
    color = :green, linestyle = :dot)
plot!(p, x_high_opt, T_opt, label = "", linewidth = 2,
    color = :green, linestyle = :dot)

# Reference data points
ref_T = [d[1] for d in REFERENCE_DATA]
ref_x_low = [d[2] for d in REFERENCE_DATA]
ref_x_high = [d[3] for d in REFERENCE_DATA]

scatter!(p, ref_x_low, ref_T, label = "Reference data", markersize = 8,
    color = :black, markershape = :circle)
scatter!(p, ref_x_high, ref_T, label = "", markersize = 8,
    color = :black, markershape = :circle)

# Save
output_path = joinpath(@__DIR__, "621_parameter_optimization.png")
savefig(p, output_path)
println("Saved to: $output_path")

# =============================================================================
# Combined Figure for JOSS Paper
# =============================================================================

println("\n" * "-" ^ 70)
println("Generating combined figure (for JOSS paper)...")
println("-" ^ 70)

# Build history by re-running with increasing iteration limits
history = []
for max_iter in [1, 2, 3, 5, 10, 20]
    res = optimize(calculate_loss, initial_params, NelderMead(),
        Optim.Options(iterations = max_iter))
    push!(history, copy(Optim.minimizer(res)))
end
push!(history, copy(Optim.minimizer(result_nm)))

# --- Panel (a): Convergence Progress ---
iter_labels = ["Iter 1", "Iter 2", "Iter 3", "Iter 5", "Iter 10", "Iter 20", "Final"]
n_hist = length(history)

pa = plot(
    xlabel = "x(Ag)",
    ylabel = "Temperature [K]",
    title = "(a) Miscibility Gap Evolution",
    legend = :bottomleft,
    xlims = (0, 1),
    ylims = (750, 1350),
    legendfontsize = 7,
    titlefontsize = 10,
)

# Color gradient: red → blue
for idx in 1:n_hist
    params = history[idx]
    fcc_hist = create_fcc_phase((params[1], params[2]), (params[3], params[4]))
    T_vals, x_low, x_high = Float64[], Float64[], Float64[]
    for T in 800.0:25.0:1300.0
        gap = find_miscibility_gap(fcc_hist, T)
        if !isnothing(gap)
            push!(T_vals, T)
            push!(x_low, gap.x1)
            push!(x_high, gap.x2)
        end
    end

    t = (idx - 1) / (n_hist - 1)
    color = RGB(1-t, 0.3*(1-t), t)
    label = iter_labels[idx]
    lw = idx == n_hist ? 2.5 : 1.5
    alpha = idx == n_hist ? 1.0 : 0.6

    plot!(pa, x_low, T_vals, label = label, linewidth = lw, color = color, alpha = alpha)
    plot!(pa, x_high, T_vals, label = "", linewidth = lw, color = color, alpha = alpha)
end

# Reference data
scatter!(
    pa,
    ref_x_low,
    ref_T,
    label = "Target",
    markersize = 7,
    color = :black,
    markershape = :star5,
)
scatter!(
    pa,
    ref_x_high,
    ref_T,
    label = "",
    markersize = 7,
    color = :black,
    markershape = :star5,
)

# --- Panel (b): Loss Convergence ---
losses = [calculate_loss(h) for h in history]
iter_nums_plot = [1, 2, 3, 5, 10, 20, Optim.iterations(result_nm)]

pb = plot(
    xlabel = "Iteration",
    ylabel = "Loss (SSE)",
    title = "(b) Optimization Convergence",
    xscale = :log10,
    yscale = :log10,
    xlims = (1, 100),
    legend = false,
    titlefontsize = 10,
)

plot!(pb, iter_nums_plot, losses, linewidth = 2, color = :gray, alpha = 0.5)

for idx in 1:n_hist
    t = (idx - 1) / (n_hist - 1)
    color = RGB(1-t, 0.3*(1-t), t)
    scatter!(
        pb,
        [iter_nums_plot[idx]],
        [losses[idx]],
        markersize = 10,
        color = color,
        markerstrokewidth = 1,
    )
end

min_loss = minimum(losses)
max_loss = maximum(losses)
improvement = round(Int, log10(max_loss / min_loss))
annotate!(
    pb,
    5,
    sqrt(min_loss * max_loss),
    text("$(improvement) orders\nimproved", 9, :center),
)

# --- Combine ---
p_combined = plot(pa, pb, layout = (1, 2), size = (900, 400), margin = 5Plots.mm)

output_combined = joinpath(@__DIR__, "621_parameter_optimization_combined.png")
savefig(p_combined, output_combined)
println("Saved to: $output_combined")

display(p_combined)

# =============================================================================
# Summary
# =============================================================================

println("\n" * "=" ^ 70)
println("Summary")
println("=" ^ 70)

println("""

This example demonstrates:

1. **Julia DSL**: Thermodynamic models defined as Julia functions
   - set_G!(fcc, [:AG], ghser_ag)
   - set_L!(fcc, [:AG, :CU], 0, T -> L0_a + L0_b*T)

2. **Parameter Optimization**: Fitting L0, L1 to phase boundary data
   - Loss function involves iterative equilibrium calculations
   - find_miscibility_gap() is called at each temperature

3. **Why this is unique to Julia + AD**:
   - Symbolic differentiation cannot differentiate through iterative solvers
   - Numerical differentiation is slow and inaccurate
   - AD provides exact gradients through the entire computation

Note: Full AD-based optimization (autodiff=:forward in Optim.jl) would
require making find_miscibility_gap differentiable, which is future work.
Currently we use numerical gradients with BFGS, but the Julia DSL
infrastructure is in place for full AD support.
""")
