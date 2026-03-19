# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of OpenCALPHAD.jl - Example: Fe-Ni Gibbs energy with magnetic contribution
#
# Plot FCC_A1 Gibbs energy curve for the Fe-Ni system,
# showing the effect of the magnetic contribution
# (Inden-Hillert-Jarl model).
#
# Uses FENI_clean.TDB (SSOL2-based, MQ& parameters removed).
#
# Usage:
#   julia --project examples/103_gibbs_energy_feni.jl

using OpenCALPHAD
using Plots

# Load Fe-Ni thermodynamic database
tdb_path = joinpath(@__DIR__, "..", "reftest", "tdb", "FENI_clean.TDB")
db = read_tdb(tdb_path)
fcc = get_phase(db, "FCC_A1")

# Temperature for calculation
T = 1000.0  # K

println("=" ^ 60)
println("Fe-Ni FCC Gibbs Energy at T=$(Int(T))K")
println("=" ^ 60)

# Scan composition to get G(x) curve
solver = GridSearchSolver(n_points = 101)
scan = scan_composition(fcc, T, db, solver)

# Calculate magnetic contribution separately
n_points = 101
x_grid = range(0.0, 1.0, length = n_points)
G_mag = Float64[]
for x in x_grid
    y = OpenCALPHAD.make_y_matrix(fcc, x)
    push!(G_mag, calculate_magnetic_energy(fcc, T, y, db))
end

println("\nMagnetic contribution at selected compositions:")
for x in [0.0, 0.3, 0.5, 0.7, 1.0]
    y = OpenCALPHAD.make_y_matrix(fcc, x)
    gm = calculate_magnetic_energy(fcc, T, y, db)
    gt = calculate_gibbs_energy(fcc, T, y, db)
    println(
        "  x(Ni)=$(x): G_mag=$(round(gm, digits=1)) J/mol, " *
        "G_total=$(round(gt, digits=1)) J/mol",
    )
end

# Plot 1: Total Gibbs energy
p1 = plot(
    scan.x_grid, scan.G_values ./ 1000,
    xlabel = "x(Ni)",
    ylabel = "Gibbs Energy [kJ/mol]",
    title = "Fe-Ni FCC at T=$(Int(T))K",
    linewidth = 2,
    label = "G_total",
    color = :blue,
    legend = :topright,
    size = (600, 400),
)

# Plot 2: Magnetic contribution
p2 = plot(
    collect(x_grid), G_mag,
    xlabel = "x(Ni)",
    ylabel = "G_magnetic [J/mol]",
    title = "Magnetic Contribution (IHJ)",
    linewidth = 2,
    label = "G_mag at $(Int(T))K",
    color = :red,
    legend = :topright,
    size = (600, 400),
)

# Combined plot
p_combined = plot(p1, p2, layout = (1, 2), size = (1100, 400),
    margin = 5Plots.mm)

output_path = joinpath(@__DIR__, "103_gibbs_energy_feni.png")
savefig(p_combined, output_path)
println("\nSaved: $output_path")

# Fortran reference check
println("\n" * "-" ^ 60)
println("Fortran reference (oc6P v6.100):")
println("  T=1000K, x(Ni)=0.3: G = -49189.449 J/mol")
y_ref = OpenCALPHAD.make_y_matrix(fcc, 0.3)
G_ref = calculate_gibbs_energy(fcc, 1000.0, y_ref, db)
println("  Julia:              G = $(round(G_ref, digits=3)) J/mol")
println("  Difference: $(round(abs(G_ref - (-49189.449)), digits=3)) J/mol")

display(p_combined)
