@testset "KL Exchange" begin
    random_seed = 443591

    Random.seed!(random_seed)
    candidates = FullFactorial(fill([-1, 0, 1], 5))
    candidates_formula = ConstantTerm(0) ~ sum(Term.(Symbol.(names(candidates.matrix))))
    candidates_formula = FormulaTerm(candidates_formula.lhs,
                                     candidates_formula.rhs +
                                     (@formula 0 ~ factor3 ^ 2).rhs)
    selected_rows = kl_exchange(candidates_formula,
                                candidates.matrix,
                                seed_design_size = 2,
                                experiments = 11,
                                design_k = 11,
                                candidates_l = size(candidates.matrix, 1) - 11)

    @test isapprox(d_criterion(selected_rows), 0.730476820204)

    Random.seed!(random_seed)
    n_candidates = 3000
    n_factors = 30;
    design_generator = DesignDistribution(Distributions.Uniform(0, 1), n_factors)
    candidates = rand(design_generator, n_candidates)
    selected_rows = kl_exchange(candidates.formula,
                                candidates.matrix,
                                experiments = n_factors + 2,
                                design_k = n_factors,
                                candidates_l = round(Int, (n_candidates - n_factors) / 2))

    @test isapprox(d_criterion(selected_rows), 0.079049188142)
end
