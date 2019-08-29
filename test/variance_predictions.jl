using Random, Logging

@testset "variance_predictions" begin
    function test_optimize_design(;factors::Int,
                                  levels::Int,
                                  experiments::Int,
                                  candidate_set_size::Int,
                                  iterations::Int,
                                  seed::Int,
                                  refresh_candidate_set::Bool,
                                  full_factorial::Bool)
        Random.seed!(seed)

        if refresh_candidate_set || !full_factorial
            design = generate_random_design(factors, experiments)
            candidate_set = generate_random_design(factors, candidate_set_size)
        else
            candidate_set = expanded_design(factors, levels)
            removed_indices = randperm(levels ^ factors)[1:experiments]

            design = candidate_set[removed_indices, :]
            candidate_set = candidate_set[[i for i in 1:(levels ^ factors) if !(i in removed_indices)], :]
        end

        optimized_design = optimize_design(factors = factors,
                                           levels = levels,
                                           experiments = experiments,
                                           design = design,
                                           candidate_set = candidate_set,
                                           iterations = iterations,
                                           refresh_candidate_set = refresh_candidate_set)
    end

    output = test_optimize_design(factors               = 4,
                                  levels                = 2,
                                  experiments           = 12,
                                  candidate_set_size    = 1000,
                                  iterations            = 40,
                                  seed                  = 403,
                                  refresh_candidate_set = false,
                                  full_factorial        = false)

    @test isapprox(d_criterion(output), 0.2605604362493956)

    output = test_optimize_design(factors               = 4,
                                  levels                = 2,
                                  experiments           = 12,
                                  candidate_set_size    = 1000,
                                  iterations            = 40,
                                  seed                  = 403,
                                  refresh_candidate_set = true,
                                  full_factorial        = false)

    @test isapprox(d_criterion(output), 0.2843243135328792)

    expected_output = [-1.0 -1.0 1.0 1.0 1.0;
                       -1.0 1.0 1.0 -1.0 1.0;
                       -1.0 -1.0 -1.0 1.0 1.0;
                       1.0 1.0 -1.0 1.0 1.0;
                       1.0 -1.0 -1.0 1.0 1.0;
                       -1.0 1.0 -1.0 1.0 1.0]

    output = test_optimize_design(factors               = 4,
                                  levels                = 2,
                                  experiments           = 6,
                                  candidate_set_size    = -1,
                                  iterations            = -1,
                                  seed                  = 403,
                                  refresh_candidate_set = false,
                                  full_factorial        = true)

    @test isapprox(d_criterion(output), 0.6666666666666666)
    @test output == expected_output
end
