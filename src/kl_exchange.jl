using LinearAlgebra, StatsModels, ExperimentalDesign, Distributions, DataFrames

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> 1 + 1
2
```
"""
function d_criterion(model_matrix;
                     tolerance = 1e-9)
    det(((model_matrix' * model_matrix) + (I * tolerance)) /
        size(model_matrix, 1)) ^
        (1 / size(model_matrix, 2))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> candidates = FullFactorial(fill([-1, 0, 1], 5));

julia> candidates_formula = ConstantTerm(0) ~ sum(Term.(names(candidates.matrix)));

julia> candidates_formula = FormulaTerm(candidates_formula.lhs,
                                        candidates_formula.rhs +
                                        (@formula 0 ~ factor3 ^ 2).rhs);

julia> selected_rows = kl_exchange(candidates_formula,
                                   candidates.matrix,
                                   seed_design_size = 2,
                                   experiments = 11,
                                   design_k = 11,
                                   candidates_l = size(candidates.matrix, 1) - 11);

julia> d_criterion(selected_rows)
0.730476820204009

julia> candidates.matrix[selected_rows.indices[1], :]
11×5 DataFrames.DataFrame
│ Row │ factor1 │ factor2 │ factor3 │ factor4 │ factor5 │
│     │ Int64   │ Int64   │ Int64   │ Int64   │ Int64   │
├─────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ 1   │ -1      │ 1       │ -1      │ 1       │ 1       │
│ 2   │ 1       │ -1      │ 1       │ 1       │ 1       │
│ 3   │ 1       │ 1       │ 0       │ -1      │ 1       │
│ 4   │ -1      │ 1       │ 1       │ -1      │ 1       │
│ 5   │ -1      │ -1      │ 0       │ 1       │ 1       │
│ 6   │ 1       │ 1       │ 1       │ 1       │ -1      │
│ 7   │ -1      │ -1      │ 1       │ -1      │ -1      │
│ 8   │ 1       │ -1      │ 0       │ 1       │ -1      │
│ 9   │ -1      │ 1       │ -1      │ 1       │ -1      │
│ 10  │ -1      │ 1       │ 0       │ -1      │ -1      │
│ 11  │ 1       │ -1      │ -1      │ -1      │ -1      │
```
"""
function kl_exchange(formula::FormulaTerm,
                     candidate_set::DataFrame;
                     tolerance = 1e-9,
                     seed_design_size = 2,
                     max_iterations::Int = 1000,
                     experiments::Int,
                     design_k::Int,
                     candidates_l::Int)
    model_matrix = ModelMatrix(ModelFrame(formula, candidate_set)).m

    seed = view(model_matrix,
                rand(1:size(model_matrix, 1), seed_design_size),
                :)

    l_variances = zeros(Int, candidates_l)
    k_variances = zeros(Int, design_k)

    design_information = zeros(size(model_matrix, 2),
                               size(model_matrix, 2))

    variances = zeros(size(model_matrix, 1) - size(seed, 1))

    # Seed Iterations
    for i = 1:(experiments - seed_design_size)
        seed_information = view(design_information, 1:size(seed, 2), 1:size(seed, 2))
        seed_information .= inv((seed' * seed) + (I * tolerance))

        candidates = view(model_matrix,
                          setdiff(1:size(model_matrix, 1), seed.indices[1]),
                          :)

        variances_candidates = view(variances, 1:size(candidates, 1))
        variances_candidates .= dot.(eachrow(candidates * seed_information),
                                     eachrow(candidates))

        l_variances .= partialsortperm(variances_candidates, 1:candidates_l, rev = true)

        seed = view(model_matrix,
                    vcat(candidates.indices[1][findmax(l_variances)[2]],
                         seed.indices[1]),
                    :)

    end

    design = seed
    variances_design = zeros(size(design, 1))
    variances_candidates = view(variances, 1:(size(model_matrix, 1) - size(design, 1)))

    kl_deltas = zeros(design_k, candidates_l)

    # Exchange Steps
    for i = 1:max_iterations
        design_information .= inv((design' * design) + (I * tolerance))
        candidates = view(model_matrix,
                          setdiff(1:size(model_matrix, 1), design.indices[1]),
                          :)

        variances_design .= dot.(eachrow(design * design_information),
                                 eachrow(design))
        k_variances .= partialsortperm(variances_design, 1:design_k, rev = false)

        variances_candidates .= dot.(eachrow(candidates * design_information),
                                     eachrow(candidates))
        l_variances .= partialsortperm(variances_candidates, 1:candidates_l, rev = true)

        k_swappable = view(model_matrix,
                           [design.indices[1][k] for k in k_variances],
                           :)

        l_swappable = view(model_matrix,
                           [candidates.indices[1][l] for l in l_variances],
                           :)

        kl_deltas .= k_swappable * design_information * l_swappable'

        swap = ((-1, -1), -1)

        for k = 1:design_k
            for l = 1:candidates_l
                delta = ((1 - variances_design[k_variances[k]]) *
                         (1 + variances_candidates[l_variances[l]])) +
                         (kl_deltas[k, l] ^ 2)

                if delta > swap[2]
                    swap = ((k, l), delta)
                end
            end
        end

        if swap[2] - 1.0 <= tolerance
            break
        end

        design = view(model_matrix,
                      vcat(deleteat!(copy(design.indices[1]),
                                     k_variances[swap[1][1]]),
                           candidates.indices[1][l_variances[swap[1][2]]]),
                      :)
    end

    return design
end

function test()
    c = 300000
    n = 50

    println("Allocating memory")
    design_generator = RandomDesign(Distributions.Uniform(0, 1), n)
    @time candidates = rand(design_generator, c)

    println("Design created, calling exchange")
    println(design_generator.formula)

    @time b = kl_exchange(design_generator.formula,
                          candidates,
                          experiments = n + 2,
                          design_k = n - 3,
                          candidates_l = round(Int, (c - n) / 2))
    return b
end

function test2()
    candidates = FullFactorial(fill([-1, 1], 10))

    println(typeof(candidates.formula.rhs))
    println(sum(Term.(names(candidates.matrix))))
    candidates_formula = ConstantTerm(0) ~ sum(Term.(names(candidates.matrix)))

    return kl_exchange(candidates_formula,
                       candidates.matrix,
                       seed_design_size = 8,
                       experiments = 11,
                       design_k = 11,
                       candidates_l = size(candidates.matrix, 1) - 11)
end

function test3()
end
