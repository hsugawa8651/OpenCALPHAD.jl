using Test
using OpenCALPHAD

@testset "OpenCALPHAD.jl" begin
    @testset "Types" begin
        @testset "Element" begin
            ag = Element("AG", "FCC_A1", 107.8682, 5744.6, 42.551)
            @test ag.name == "AG"
            @test ag.reference_state == "FCC_A1"
            @test ag.mass ≈ 107.8682
        end

        @testset "Phase" begin
            fcc = Phase("FCC_A1", 'S', [1.0])
            @test fcc.name == "FCC_A1"
            @test fcc.state == 'S'
            @test num_sublattices(fcc) == 1
            @test total_sites(fcc) ≈ 1.0
        end

        @testset "Conditions" begin
            conds = ConditionSet()
            @test conds.temperature ≈ 298.15
            @test conds.pressure ≈ 1e5

            set_condition!(conds, EquilCondition(COND_T, 1000.0))
            @test conds.temperature ≈ 1000.0

            set_condition!(conds, EquilCondition(COND_X, "CU", 0.3))
            @test conds.compositions["CU"] ≈ 0.3
        end
    end

    @testset "Database" begin
        db = Database()
        @test length(db.elements) == 0
        @test length(db.phases) == 0
    end

    @testset "Constants" begin
        @test OpenCALPHAD.R ≈ 8.314462618
        @test OpenCALPHAD.P_REF ≈ 1e5
    end

    @testset "Minimizer" begin
        using LinearAlgebra: norm

        # Setup: Create a simple binary phase for testing
        function create_test_phase_and_db()
            db = Database()

            # Add elements
            push!(db.elements, Element("CR", "BCC_A2", 51.996, 4050.0, 23.543))
            push!(db.elements, Element("FE", "BCC_A2", 55.845, 4489.0, 27.28))

            # Create BCC phase with single sublattice (CR,FE)
            bcc = Phase("BCC_A2", 'S', [1.0])
            push!(bcc.constituents, ["CR", "FE"])
            push!(db.phases, bcc)

            # Add simple Gibbs energy parameters (simplified for testing)
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

        @testset "Fallback errors" begin
            bcc, db = create_test_phase_and_db()

            # Test unsupported solver type
            struct FakeSolver end
            @test_throws ErrorException minimize_gibbs(bcc, 1000.0, db, FakeSolver())

            # Test unsupported diff method
            struct FakeDiffMethod end
            y = [0.5 0.5]
            @test_throws ErrorException compute_gradient(
                bcc,
                1000.0,
                y,
                db,
                FakeDiffMethod(),
            )
            @test_throws ErrorException compute_hessian(
                bcc,
                1000.0,
                y,
                db,
                FakeDiffMethod(),
            )

            # Test invalid phase type for calculate_gibbs_energy
            @test_throws ErrorException calculate_gibbs_energy("not a phase", 1000.0, y, db)
            @test_throws ErrorException calculate_gibbs_energy("not a phase", 1000.0, y)
        end

        @testset "GridSearchSolver" begin
            bcc, db = create_test_phase_and_db()
            solver = GridSearchSolver(n_points = 11)

            result = scan_composition(bcc, 1000.0, db, solver)

            @test result isa GridScanResult
            @test result.phase === bcc
            @test length(result.x_grid) == 11
            @test length(result.G_values) == 11
            @test length(result.y_values) == 11
            @test 1 <= result.min_idx <= 11

            # G values should be finite
            @test all(isfinite, result.G_values)
        end

        @testset "Differentiation methods" begin
            bcc, db = create_test_phase_and_db()
            y = [0.3 0.7]  # 30% CR, 70% FE

            # Test ForwardDiff gradient
            grad_fd = compute_gradient(bcc, 1000.0, y, db, ForwardDiffMethod())
            @test length(grad_fd) == 2
            @test all(isfinite, grad_fd)

            # Test FiniteDiff gradient
            grad_fin = compute_gradient(bcc, 1000.0, y, db, FiniteDiffMethod())
            @test length(grad_fin) == 2
            @test all(isfinite, grad_fin)

            # Compare: should be close (relative error < 1e-5)
            rel_error = norm(grad_fd - grad_fin) / max(norm(grad_fd), 1.0)
            @test rel_error < 1e-5

            # Test Hessians
            hess_fd = compute_hessian(bcc, 1000.0, y, db, ForwardDiffMethod())
            hess_fin = compute_hessian(bcc, 1000.0, y, db, FiniteDiffMethod())

            @test size(hess_fd) == (2, 2)
            @test size(hess_fin) == (2, 2)

            # Hessian comparison: just check both are finite and symmetric
            # (finite diff Hessian can be inaccurate with log terms)
            @test all(isfinite, hess_fd)
            @test all(isfinite, hess_fin)
            @test hess_fd ≈ hess_fd'  # Symmetric
            @test hess_fin ≈ hess_fin'  # Symmetric
        end

        @testset "Chemical potential" begin
            bcc, db = create_test_phase_and_db()
            y = [0.5 0.5]

            mu = compute_chemical_potential(bcc, 1000.0, y, db, ForwardDiffMethod())
            @test length(mu) == 2
            @test all(isfinite, mu)

            # At equimolar composition, chemical potentials should differ
            # (unless it's an ideal solution at special T)
            # Just check they are reasonable values (order of 1e4 to 1e5 J/mol)
            @test all(abs.(mu) .< 1e7)
        end

        @testset "NewtonSolver" begin
            bcc, db = create_test_phase_and_db()

            # Create Newton solver
            newton = NewtonSolver(
                diff_method = ForwardDiffMethod(),
                max_iter = 100,
                ftol = 1e-8,
            )

            # Test minimize_gibbs returns SinglePhaseResult
            result = minimize_gibbs(bcc, 1000.0, db, newton)

            @test result isa SinglePhaseResult
            @test result.phase === bcc
            @test size(result.site_fractions) == (1, 2)
            @test isfinite(result.gibbs_energy)
            @test result.converged == true
            @test result.iterations > 0

            # Site fractions should sum to 1
            @test sum(result.site_fractions[1, :]) ≈ 1.0 atol=1e-10

            # Site fractions should be in [0, 1]
            @test all(0 .<= result.site_fractions .<= 1)
        end

        @testset "HybridSolver" begin
            bcc, db = create_test_phase_and_db()

            # Create Hybrid solver (grid + Newton)
            hybrid = HybridSolver(
                n_points = 11,
                diff_method = ForwardDiffMethod(),
                max_iter = 100,
                ftol = 1e-8,
            )

            result = minimize_gibbs(bcc, 1000.0, db, hybrid)

            @test result isa SinglePhaseResult
            @test result.converged == true

            # Hybrid should give same or better result than grid alone
            grid_solver = GridSearchSolver(n_points = 11)
            grid_result = scan_composition(bcc, 1000.0, db, grid_solver)
            grid_min_G = grid_result.G_values[grid_result.min_idx]

            @test result.gibbs_energy <= grid_min_G + 1.0  # Newton should refine
        end

        @testset "clamp_site_fractions!" begin
            y = [0.0 1.0; 1e-15 0.9999999999999985]

            clamp_site_fractions!(y; min_val = 1e-10)

            # Check minimum value constraint
            @test all(y .>= 1e-10)

            # Check normalization (each row sums to 1)
            @test sum(y[1, :]) ≈ 1.0
            @test sum(y[2, :]) ≈ 1.0
        end
    end

    @testset "Two-phase equilibrium" begin
        using LinearAlgebra: norm

        # Setup: Create two phases (BCC and FCC) for Cr-Fe
        function create_two_phase_system()
            db = Database()

            # Add elements
            push!(db.elements, Element("CR", "BCC_A2", 51.996, 4050.0, 23.543))
            push!(db.elements, Element("FE", "BCC_A2", 55.845, 4489.0, 27.28))

            # Create BCC phase
            bcc = Phase("BCC_A2", 'S', [1.0])
            push!(bcc.constituents, ["CR", "FE"])
            push!(db.phases, bcc)

            # BCC parameters
            push!(
                bcc.parameters,
                Parameter(
                    "BCC_A2", 'G', [["CR"]], 0,
                    [298.15, 6000.0],
                    ["-8856.94+157.48*T-26.908*T*LN(T)"],
                ),
            )
            push!(
                bcc.parameters,
                Parameter(
                    "BCC_A2", 'G', [["FE"]], 0,
                    [298.15, 6000.0],
                    ["1225.7+124.134*T-23.5143*T*LN(T)"],
                ),
            )
            push!(
                bcc.parameters,
                Parameter(
                    "BCC_A2", 'L', [["CR", "FE"]], 0,
                    [298.15, 6000.0],
                    ["20500-9.68*T"],
                ),
            )

            # Create FCC phase
            fcc = Phase("FCC_A1", 'S', [1.0])
            push!(fcc.constituents, ["CR", "FE"])
            push!(db.phases, fcc)

            # FCC parameters (shifted to create a two-phase region)
            push!(
                fcc.parameters,
                Parameter(
                    "FCC_A1", 'G', [["CR"]], 0,
                    [298.15, 6000.0],
                    ["-1572.94+157.48*T-26.908*T*LN(T)"],  # Higher than BCC
                ),
            )
            push!(
                fcc.parameters,
                Parameter(
                    "FCC_A1", 'G', [["FE"]], 0,
                    [298.15, 6000.0],
                    ["-1462.4+124.134*T-23.5143*T*LN(T)"],  # Lower than BCC
                ),
            )
            push!(
                fcc.parameters,
                Parameter(
                    "FCC_A1", 'L', [["CR", "FE"]], 0,
                    [298.15, 6000.0],
                    ["10833-7.477*T"],  # Different interaction
                ),
            )

            return bcc, fcc, db
        end

        @testset "TwoPhaseResult type" begin
            bcc, fcc, db = create_two_phase_system()

            # Manually create a TwoPhaseResult for type testing
            y1 = [0.3 0.7]
            y2 = [0.7 0.3]
            G1 = calculate_gibbs_energy(bcc, 1000.0, y1, db)
            G2 = calculate_gibbs_energy(fcc, 1000.0, y2, db)

            r1 = SinglePhaseResult(bcc, y1, G1, true, 10)
            r2 = SinglePhaseResult(fcc, y2, G2, true, 10)

            result = TwoPhaseResult(r1, r2, 0.4, 0.6, [G1, G2], true, 10)

            @test result.phase1 === r1
            @test result.phase2 === r2
            @test result.phase1_fraction ≈ 0.4
            @test result.phase2_fraction ≈ 0.6
            @test result.converged == true
        end
    end

    @testset "STEP calculation" begin
        # Load Ag-Cu database for testing
        tdb_path = joinpath(dirname(@__DIR__), "reftest", "tdb", "agcu.TDB")
        db = read_tdb(tdb_path)
        fcc = get_phase(db, "FCC_A1")

        @testset "Type error fallbacks" begin
            # Test invalid phase type
            @test_throws ErrorException step_temperature(
                "not a phase", db, 0.7, 800.0, 1000.0, 100.0,
            )

            # Test invalid database type
            @test_throws ErrorException step_temperature(
                fcc, "not a database", 0.7, 800.0, 1000.0, 100.0,
            )
        end

        @testset "StepPoint construction" begin
            # Create a StepPoint manually
            point = StepPoint(
                1000.0,           # axis_value (T)
                -53648.0,         # gibbs_energy
                2,                # n_phases
                [0.228, 0.772],   # phase_amounts
                [[0.034, 0.966], [0.897, 0.103]],  # phase_compositions
                true,              # converged
            )

            @test point.axis_value ≈ 1000.0
            @test point.gibbs_energy ≈ -53648.0
            @test point.n_phases == 2
            @test length(point.phase_amounts) == 2
            @test point.phase_amounts[1] ≈ 0.228
            @test point.phase_amounts[2] ≈ 0.772
            @test point.converged == true
        end

        @testset "StepResult construction" begin
            # Create StepPoints
            point1 = StepPoint(
                800.0,
                -39427.0,
                2,
                [0.275, 0.725],
                [[0.009, 0.991], [0.962, 0.038]],
                true,
            )
            point2 = StepPoint(
                900.0,
                -46335.0,
                2,
                [0.256, 0.744],
                [[0.018, 0.982], [0.935, 0.065]],
                true,
            )
            point3 = StepPoint(
                1000.0,
                -53648.0,
                2,
                [0.228, 0.772],
                [[0.034, 0.966], [0.897, 0.103]],
                true,
            )

            # Create StepResult
            result = StepResult(:T, [800.0, 900.0, 1000.0], [point1, point2, point3], 0.7)

            @test result.axis_variable == :T
            @test length(result.axis_values) == 3
            @test length(result.points) == 3
            @test result.x_overall ≈ 0.7
            @test result.points[1].axis_value ≈ 800.0
            @test result.points[3].gibbs_energy ≈ -53648.0
        end
    end

    @testset "MAP calculation (phase diagram)" begin
        # Load Ag-Cu database for testing
        tdb_path = joinpath(dirname(@__DIR__), "reftest", "tdb", "agcu.TDB")
        db = read_tdb(tdb_path)
        fcc = get_phase(db, "FCC_A1")

        @testset "Type error fallbacks" begin
            # Test invalid phase type
            @test_throws ErrorException map_phase_diagram(
                "not a phase", db, 800.0, 1000.0, 100.0,
            )

            # Test invalid database type
            @test_throws ErrorException map_phase_diagram(
                fcc, "not a database", 800.0, 1000.0, 100.0,
            )
        end

        @testset "PhaseDiagramPoint construction" begin
            # Create a PhaseDiagramPoint manually
            point = PhaseDiagramPoint(
                1000.0,           # temperature
                (0.034, 0.897),   # compositions tuple (Cu-rich, Ag-rich)
                true,              # converged
            )

            @test point.temperature ≈ 1000.0
            @test point.compositions[1] ≈ 0.034
            @test point.compositions[2] ≈ 0.897
            @test point.converged == true
        end

        @testset "PhaseDiagramResult construction" begin
            # Create points
            point1 = PhaseDiagramPoint(800.0, (0.009, 0.962), true)
            point2 = PhaseDiagramPoint(900.0, (0.018, 0.935), true)
            point3 = PhaseDiagramPoint(1000.0, (0.034, 0.897), true)

            # Create result
            result = PhaseDiagramResult(
                fcc,
                [point1, point2, point3],
            )

            @test result.phase === fcc
            @test length(result.points) == 3
            @test result.points[1].temperature ≈ 800.0
            @test result.points[3].compositions[2] ≈ 0.897
        end
    end

    # Multi-component support tests (Cr-Mo from steel1.TDB)
    include("multicomponent.jl")
end
