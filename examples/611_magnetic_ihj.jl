# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Inden-Hillert-Jarl magnetic model
#
# Inden-Hillert-Jarl Magnetic Contribution Example
#
# Demonstrates conditional branching in Gibbs energy functions
# using the magnetic contribution model for BCC Fe.
#
# This shows a use case where Julia's if/else is more natural
# than symbolic Piecewise expressions.
#
# References:
#   [1] M. Hillert and M. Jarl, "A model for alloying in ferromagnetic metals",
#       Calphad, vol. 2, no. 3, pp. 227-238, 1978.
#       https://doi.org/10.1016/0364-5916(78)90011-1
#
#   [2] G. Inden, "Project Meeting CALPHAD V", Max-Planck-Inst.
#       Eisenforschung, Düsseldorf, Germany, Chapter IX, 1976.
#
#   [3] A.T. Dinsdale, "SGTE data for pure elements",
#       Calphad, vol. 15, no. 4, pp. 317-425, 1991.
#       https://doi.org/10.1016/0364-5916(91)90030-N
#
# Usage:
#   julia --project examples/611_magnetic_ihj.jl

using OpenCALPHAD
using Plots
using ForwardDiff

println("=" ^ 70)
println("Inden-Hillert-Jarl Magnetic Contribution Model")
println("=" ^ 70)

# =============================================================================
# Physical constants and parameters
# =============================================================================

const R = 8.314462618  # Gas constant [J/(mol·K)]

# BCC Fe parameters
const Tc_Fe = 1043.0   # Curie temperature [K]
const beta_Fe = 2.22   # Magnetic moment [Bohr magneton]
const p_bcc = 0.40     # Structure factor for BCC

# =============================================================================
# Inden-Hillert-Jarl Model
# =============================================================================

"""
    inden_hillert_g(tau, p)

Compute the g(τ) function in the Inden-Hillert-Jarl model.

  - τ < 1: Ferromagnetic region (below Curie temperature)
  - τ ≥ 1: Paramagnetic region (above Curie temperature)

The function has different expressions in each region,
demonstrating natural use of conditional branching.
"""
function inden_hillert_g(tau, p)
    # Normalization constant
    A = 518/1125 + (11692/15975) * (1/p - 1)

    if tau < 1
        # Ferromagnetic region (T < Tc)
        # Polynomial expansion for τ < 1
        term1 = 79 / (140 * p)
        term2 = (474/497) * (1/p - 1) * (tau^3/6 + tau^9/135 + tau^15/600)
        g = 1 - (term1 * tau^(-1) + term2) / A
    else
        # Paramagnetic region (T > Tc)
        # Inverse power series for τ ≥ 1
        g = -(tau^(-5)/10 + tau^(-15)/315 + tau^(-25)/1500) / A
    end

    return g
end

"""
    G_magnetic(T, Tc, beta, p)

Compute the magnetic contribution to Gibbs energy [J/mol].

G_mag = RT * ln(β + 1) * g(τ)

where τ = T/Tc is the reduced temperature.
"""
function G_magnetic(T, Tc, beta, p)
    tau = T / Tc
    g = inden_hillert_g(tau, p)
    return R * T * log(beta + 1) * g
end

# =============================================================================
# GHSER for BCC Fe (without magnetic contribution)
# =============================================================================

function ghser_fe_nonmag(T)
    # Simplified GHSER for Fe (non-magnetic part)
    # From SGTE database
    if T < 1811.0  # Below melting point
        return 1225.7 + 124.134*T - 23.5143*T*log(T) -
               0.00439752*T^2 - 5.89269e-8*T^3 + 77358.5/T
    else
        return -25383.581 + 299.31255*T - 46.0*T*log(T) +
               2.2960305e31*T^(-9)
    end
end

# =============================================================================
# Full GHSER for BCC Fe (with magnetic contribution)
# =============================================================================

function ghser_fe_full(T)
    return ghser_fe_nonmag(T) + G_magnetic(T, Tc_Fe, beta_Fe, p_bcc)
end

# =============================================================================
# Verify automatic differentiation works through conditional
# =============================================================================

println("\nVerifying ForwardDiff through conditional branching...")

# Test at T < Tc (ferromagnetic)
T_ferro = 800.0
G_ferro = G_magnetic(T_ferro, Tc_Fe, beta_Fe, p_bcc)
dG_ferro = ForwardDiff.derivative(T -> G_magnetic(T, Tc_Fe, beta_Fe, p_bcc), T_ferro)
println("  T = $T_ferro K (τ = $(T_ferro/Tc_Fe), ferromagnetic)")
println("    G_mag = $(round(G_ferro, digits=2)) J/mol")
println("    dG_mag/dT = $(round(dG_ferro, digits=4)) J/(mol·K)")

# Test at T > Tc (paramagnetic)
T_para = 1200.0
G_para = G_magnetic(T_para, Tc_Fe, beta_Fe, p_bcc)
dG_para = ForwardDiff.derivative(T -> G_magnetic(T, Tc_Fe, beta_Fe, p_bcc), T_para)
println("  T = $T_para K (τ = $(T_para/Tc_Fe), paramagnetic)")
println("    G_mag = $(round(G_para, digits=2)) J/mol")
println("    dG_mag/dT = $(round(dG_para, digits=4)) J/(mol·K)")

# Second derivative (heat capacity contribution)
d2G_ferro = ForwardDiff.derivative(
    T -> ForwardDiff.derivative(t -> G_magnetic(t, Tc_Fe, beta_Fe, p_bcc), T),
    T_ferro,
)
println("  d²G_mag/dT² at 800 K = $(round(d2G_ferro, digits=6)) J/(mol·K²)")

# =============================================================================
# Plot G_magnetic vs Temperature
# =============================================================================

println("\nGenerating plots...")

T_range = 300.0:10.0:1800.0
G_mag_values = [G_magnetic(T, Tc_Fe, beta_Fe, p_bcc) for T in T_range]
G_full_values = [ghser_fe_full(T) for T in T_range]
G_nonmag_values = [ghser_fe_nonmag(T) for T in T_range]

# Plot 1: Magnetic contribution only
p1 = plot(T_range, G_mag_values,
    xlabel = "Temperature [K]",
    ylabel = "G_magnetic [J/mol]",
    title = "Magnetic Contribution to Gibbs Energy (BCC Fe)",
    linewidth = 2,
    color = :blue,
    legend = false,
    size = (600, 400),
)
vline!(p1, [Tc_Fe], linestyle = :dash, color = :red, label = "")
annotate!(p1, Tc_Fe + 50, minimum(G_mag_values) * 0.5,
    text("Tc = $Tc_Fe K", 9, :left))

# Plot 2: Full vs non-magnetic GHSER
p2 = plot(T_range, G_full_values / 1000,
    label = "With magnetic",
    xlabel = "Temperature [K]",
    ylabel = "G [kJ/mol]",
    title = "GHSER for BCC Fe",
    linewidth = 2,
    color = :blue,
    size = (600, 400),
)
plot!(p2, T_range, G_nonmag_values / 1000,
    label = "Without magnetic",
    linewidth = 2,
    color = :orange,
    linestyle = :dash,
)
vline!(p2, [Tc_Fe], linestyle = :dot, color = :gray, label = "Tc")

# Plot 3: g(τ) function showing the branch
tau_range = 0.1:0.01:2.0
g_values = [inden_hillert_g(tau, p_bcc) for tau in tau_range]

p3 = plot(tau_range, g_values,
    xlabel = "τ = T/Tc",
    ylabel = "g(τ)",
    title = "Inden-Hillert g(τ) Function",
    linewidth = 2,
    color = :green,
    legend = false,
    size = (600, 400),
)
vline!(p3, [1.0], linestyle = :dash, color = :red, label = "")
annotate!(p3, 0.5, 0.5, text("Ferromagnetic\n(τ < 1)", 9, :center))
annotate!(p3, 1.5, -0.05, text("Paramagnetic\n(τ > 1)", 9, :center))

# Combine plots
p_combined = plot(p1, p3, layout = (1, 2), size = (1000, 400), margin = 5Plots.mm)

# Save
output_path = joinpath(@__DIR__, "611_magnetic_ihj.png")
savefig(p_combined, output_path)
println("Saved to: $output_path")

# =============================================================================
# Demonstrate DSL usage with magnetic contribution
# =============================================================================

println("\n" * "-" ^ 70)
println("Using Julia DSL with magnetic contribution")
println("-" ^ 70)

# Create a phase using the DSL
bcc_fe = Phase("BCC_A2", [1.0], [[:FE, :VA]])

# Set G using our function with magnetic contribution
set_G!(bcc_fe, [:FE], ghser_fe_full)
set_G!(bcc_fe, [:VA], T -> 0.0)  # Vacancy reference

println("\nPhase defined with magnetic contribution:")
println("  Phase: $(bcc_fe.name)")
println("  Julia parameters: $(length(bcc_fe.julia_parameters))")

# Verify the DSL version works
y_test = [1.0 0.0; 1.0 0.0]  # Pure Fe
G_dsl = calculate_gibbs_energy(bcc_fe, 1000.0, y_test)
println("  G at 1000 K (DSL): $(round(G_dsl, digits=2)) J/mol")

# Compare with direct calculation
G_direct = ghser_fe_full(1000.0)
println("  G at 1000 K (direct): $(round(G_direct, digits=2)) J/mol")

# =============================================================================
# Key advantage: AD through conditional branching
# =============================================================================

println("\n" * "-" ^ 70)
println("Key Advantage: Automatic Differentiation through Conditionals")
println("-" ^ 70)

println("""
In pycalphad, the magnetic model requires Piecewise symbolic expressions:

    from sympy import Piecewise, Symbol
    tau = Symbol('tau')
    g = Piecewise(
        (ferromagnetic_expr, tau < 1),
        (paramagnetic_expr, tau >= 1)
    )

In OpenCALPHAD.jl, we use natural Julia code:

    function inden_hillert_g(tau, p)
        if tau < 1
            # ferromagnetic
            ...
        else
            # paramagnetic
            ...
        end
    end

ForwardDiff.jl differentiates through both branches automatically.
This enables gradient-based optimization of magnetic parameters (Tc, β).
""")

display(p_combined)
