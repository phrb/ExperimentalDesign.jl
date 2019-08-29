using Logging, LinearAlgebra, DataFrames, StatsModels, Random

"""
$(TYPEDSIGNATURES)
"""
function d_criterion(model_matrix::Array{Float64, 2})
    det((model_matrix' * model_matrix) / size(model_matrix, 1)) ^ (1 / size(model_matrix, 2))
end

"""
$(TYPEDSIGNATURES)
"""
function linear_model(data::DataFrame)
    Term(:_) ~ sum(term.(names(data)))
end

"""
$(TYPEDSIGNATURES)

Computes `deltas` for  the exchanges between all pairs  made between experiments
in `model_matrix`  and experiments in `candidate_set`.   The `dispersion_matrix`
array   must   contain  the   dispersion   matrix   for  `model_matrix`.    Both
`model_matrix`  and   `candidate_set`  must  be  model   matrices.   The  arrays
`design_variance`,  `candidate_variance`,  `predictions`  and `deltas`  will  be
overwritten.
"""
function variance_prediction!(dispersion_matrix::Array{Float64, 2},
                              model_matrix::Array{Float64, 2},
                              candidate_set::Array{Float64, 2},
                              design_variance::Array{Float64, 1},
                              candidate_variance::Array{Float64, 1},
                              predictions::Array{Float64, 2},
                              deltas::Array{Float64, 2})
    @inbounds Threads.@threads for i = 1:size(model_matrix, 1)
        design_variance[i] = model_matrix[i, :]' * dispersion_matrix * model_matrix[i, :]
    end

    @inbounds Threads.@threads for i = 1:size(candidate_set, 1)
        candidate_variance[i] = candidate_set[i, :]' * dispersion_matrix * candidate_set[i, :]
    end

    @inbounds Threads.@threads for i = 1:size(predictions, 1)
        @inbounds for j = 1:size(predictions, 2)
            predictions[i, j] = model_matrix[i, :]' * dispersion_matrix * candidate_set[j, :]
        end
    end

    @inbounds Threads.@threads for i = 1:size(deltas, 1)
        @inbounds for j = 1:size(deltas, 2)
            deltas[i, j] = candidate_variance[j] - design_variance[i] -
                           ((design_variance[i] * candidate_variance[j]) -
                            (predictions[i, j] ^ 2))
        end
    end

    return
end

"""
$(TYPEDSIGNATURES)

Generates a random design of size `size` with `factors` numerical factors.
"""
function generate_random_design(factors::Int, size::Int)
    data = Dict("_" => ones(size))

    for i in 1:factors
        data["X$i"] = rand(size)
    end

    DataFrame(data)
end

"""
$(TYPEDSIGNATURES)

Generates  a design  from all  combinations  of `factors`  factors and  `levels`
levels.
"""
function expanded_design(factors::Int, levels::Int)
    expanded_grid = expand_grid(factors, levels)
    data = Dict("_" => ones(levels ^ factors))

    for i in 1:factors
        data["X$i"] = expanded_grid[:, i]
    end

    DataFrame(data)
end

"""
$(TYPEDSIGNATURES)

Generates all combinations of `factors` factors with `levels` levels.
"""
function expand_grid(factors::Int, levels::Int)
    chosen_levels = (((collect(1:levels) .- 1) / (levels - 1)) .* (2)) .- 1
    complete_set  = Array{Float64, 2}(undef, levels ^ factors, factors)
    iterator      = Base.product(repeat([chosen_levels], inner = factors)...)
    set_index     = 1

    for i in iterator
        complete_set[set_index, :] = collect(i)
        set_index += 1
    end

    complete_set
end

"""
$(TYPEDSIGNATURES)
"""
function optimize_design(;factors::Int,
                         levels::Int,
                         experiments::Int,
                         design::DataFrame,
                         candidate_set::DataFrame,
                         iterations::Int,
                         refresh_candidate_set::Bool = true)
    candidate_set_size = size(candidate_set, 1)

    formula = linear_model(design)

    new_schema = schema(formula, design)
    schema_formula = apply_schema(formula, new_schema)

    dummy_response, model_matrix = modelcols(schema_formula, design)
    dummy_reponse, candidate_model_matrix = modelcols(schema_formula, candidate_set)

    design_variance = similar(model_matrix[:, 1])
    candidate_variance = similar(candidate_model_matrix[:, 1])

    predictions = similar(candidate_model_matrix[:, 1],
                          (size(model_matrix, 1),
                           size(candidate_model_matrix, 1)))

    deltas = similar(candidate_model_matrix[:, 1],
                     (size(model_matrix, 1),
                      size(candidate_model_matrix, 1)))

    dispersion_matrix = inv(model_matrix' * model_matrix)

    variance_prediction!(dispersion_matrix,
                         model_matrix,
                         candidate_model_matrix,
                         design_variance,
                         candidate_variance,
                         predictions,
                         deltas)

    current_d_criterion = d_criterion(model_matrix)
    @debug "Starting D-Criterion" current_d_criterion

    for i in 1:iterations
        best_swap = findmax(deltas)
        @debug "Best Swap" best_swap

        if best_swap[1] <= 0.0
            @debug "Max. delta was less than, or equal to, zero"

            if refresh_candidate_set
                @debug "Regenerating candidate set..." i

                candidate_set = generate_random_design(factors, candidate_set_size)

                formula = linear_model(candidate_set)

                new_schema = schema(formula, candidate_set)
                schema_formula = apply_schema(formula, new_schema)

                dummy_reponse, candidate_model_matrix = modelcols(schema_formula, candidate_set)
            else
                @debug "Stopping..."
                break
            end
        else
            model_matrix[best_swap[2][1], :] = candidate_model_matrix[best_swap[2][2], :]
        end

        dispersion_matrix = inv(model_matrix' * model_matrix)

        variance_prediction!(dispersion_matrix,
                             model_matrix,
                             candidate_model_matrix,
                             design_variance,
                             candidate_variance,
                             predictions,
                             deltas)
    end

    current_d_criterion = d_criterion(model_matrix)
    @debug "Final D-Criterion" current_d_criterion

    model_matrix
end

"""
$(TYPEDSIGNATURES)
"""
function measure_optimize_design_refresh(refresh_candidate_set)
    factors = 4
    levels = 2
    experiments = 12
    candidate_set_size = 1000
    iterations = 40

    design = generate_random_design(factors, experiments)
    candidate_set = generate_random_design(factors, candidate_set_size)

    optimized_design = optimize_design(factors = factors,
                                       levels = levels,
                                       experiments = experiments,
                                       design = design,
                                       candidate_set = candidate_set,
                                       iterations = iterations,
                                       refresh_candidate_set = refresh_candidate_set)
end

"""
$(TYPEDSIGNATURES)
"""
function measure_optimize_design(levels)
    factors = 4
    experiments = 6
    iterations = 40
    refresh_candidate_set = false

    candidate_set = expanded_design(factors, levels)
    removed_indices = randperm(levels ^ factors)[1:experiments]

    design = candidate_set[removed_indices, :]
    candidate_set = candidate_set[[i for i in 1:(levels ^ factors) if !(i in removed_indices)], :]

    optimized_design = optimize_design(factors = factors,
                                       levels = levels,
                                       experiments = experiments,
                                       design = design,
                                       candidate_set = candidate_set,
                                       iterations = iterations,
                                       refresh_candidate_set = refresh_candidate_set)
end

seed = 403
Random.seed!(seed)

design = measure_optimize_design_refresh(false)
design_refresh = measure_optimize_design_refresh(true)
design_expanded = measure_optimize_design(2)
