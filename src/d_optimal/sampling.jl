function sample_full_factorial(factors::Array{T, 1}) where T <: Any
    return Array{Any, 1}([rand(i) for i in factors])
end

function check_repeated_row(subset::Array{T, 2}, row::Array{T, 1}) where T <: Any
    for j = 1:size(subset, 1)
        if subset[j, :] == row
            return true
        end
    end

    return false
end

function full_factorial_subset(factors::OrderedDict, experiments::Int)
    subset = fill!(Array{Any, 2}(experiments, length(keys(factors))), 0.0)

    for i = 1:experiments
        sample_row = sample_full_factorial(collect(values(factors)))

        while check_repeated_row(subset, sample_row)
            sample_row = sample_full_factorial(collect(values(factors)))
        end

        subset[i, :] = sample_row
    end

    subset = DataFrame(subset)
    rename!(subset, zip(names(subset), collect(keys(factors))))

    for factor in keys(factors)
        if eltype(subset[factor]) != eltype(factors[factor])
            subset[factor] = convert(Array{eltype(factors[factor]), 1}, subset[factor])
        end
    end

    return subset
end

function enforce_bounds(factors::OrderedDict,
                        sample_range::UnitRange{Int},
                        designs::Int,
                        check_bounds::Bool)
    println("> Factors: ", factors)

    full_factorial_size    = prod(length, values(factors))
    full_factorial_subsets = 2.0 ^ full_factorial_size

    println("> Full Factorial Size: ", full_factorial_size)
    println("> Total Subsets: ", full_factorial_subsets)
    println("> Range of Design Sizes: ", sample_range)
    println("> Number of Designs to Sample: ", designs)

    if check_bounds
        if sample_range.start == sample_range.stop
            restricted_subsets = factorial(float(full_factorial_size)) /
                                 (factorial(float(full_factorial_size - sample_range.start)) *
                                 factorial(float(sample_range.start)))
            println("> Total Subsets for Fixed Size ",
                    sample_range.start, ": ",
                    restricted_subsets)

            if designs > restricted_subsets
                println("> Requested too many designs, using ",
                        restricted_subsets, " instead")
                designs = floor(Int, restricted_subsets)
            end
        elseif designs > full_factorial_subsets
            println("> Requested too many designs, using ",
                    full_factorial_subsets, " instead")
            designs = floor(Int, full_factorial_subsets)
        end

        if sample_range.stop > full_factorial_size
            println("> Requested too many maximum experiments, using ",
                    full_factorial_size, " instead")
            sample_range = sample_range.start:floor(Int, full_factorial_size)
        end
    else
        println("> WARNING: Skipping bounds check!")
    end

    return designs, sample_range
end

function evaluate_all_metrics(factors::OrderedDict,
                              formula::Formula,
                              sample_range::UnitRange{Int},
                              designs::Int,
                              scale::Function)
    evaluation = DataFrame(Length  = [],
                           D       = [],
                           DELB    = [],
                           DELB_ad = [],
                           A       = [],
                           V       = [],
                           G       = [],
                           CN      = [],
                           GE      = [],
                           log2CN  = [],
                           log10D  = [])

    for i in 1:designs
        samples   = rand(sample_range)
        subset    = full_factorial_subset(factors, samples)
        candidate = generate_model_matrix(formula, subset, factors,
                                          scale = scale)

        d_opt = d_optimality(candidate)
        c_n   = condition_number(candidate)

        push!(evaluation, [size(candidate, 1),
                           d_opt,
                           d_efficiency_lower_bound(candidate),
                           d_efficiency_lower_bound_algdesign(candidate),
                           a_optimality(candidate),
                           v_optimality(candidate),
                           g_optimality(candidate),
                           c_n,
                           g_efficiency(candidate),
                           log(2, abs(c_n)),
                           log(10, abs(d_opt))])
    end

    return evaluation
end

function evaluate_metrics(factors::OrderedDict,
                          formula::Formula,
                          sample_range::UnitRange{Int},
                          designs::Int,
                          scale::Function)
    evaluation = DataFrame(Length  = [],
                           D       = [],
                           DELB    = [])

    for i in 1:designs
        samples   = rand(sample_range)
        subset    = full_factorial_subset(factors, samples)
        candidate = generate_model_matrix(formula, subset, factors,
                                          scale = scale)

        push!(evaluation, [size(candidate, 1),
                           d_optimality(candidate),
                           d_efficiency_lower_bound(candidate)])
    end

    return evaluation
end

function generate_designs(factors::OrderedDict,
                          formula::Formula,
                          sample_range::UnitRange{Int},
                          designs::Int;
                          check_bounds::Bool = true,
                          scale::Function = scale_boxdraper_encoding!,
                          compute_all_metrics::Bool = false)
    designs, sample_range = enforce_bounds(factors, sample_range,
                                           designs, check_bounds)

    if compute_all_metrics
        return evaluate_all_metrics(factors::OrderedDict,
                                    formula::Formula,
                                    sample_range::UnitRange{Int},
                                    designs::Int,
                                    scale::Function)
    else
        return evaluate_metrics(factors::OrderedDict,
                                formula::Formula,
                                sample_range::UnitRange{Int},
                                designs::Int,
                                scale::Function)
    end
end


function sample_subset(factors::OrderedDict,
                       sample_range::UnitRange{Int},
                       designs::Int;
                       check_bounds::Bool = true,
                       scale::Function = scale_boxdraper_encoding!)
    formula = build_linear_formula(collect(keys(factors)))

    # Any formula could potentially be used here:
    #   formula = @formula(y ~ x1 * x2 + x3)

    run_time = @elapsed sampling_subset = generate_designs(factors,
                                                           formula,
                                                           sample_range,
                                                           designs,
                                                           check_bounds = check_bounds,
                                                           scale = scale)
    println("> Elapsed Time: ", run_time, " seconds")

    sort!(sampling_subset, :D, rev = true)

    return sampling_subset
end

function sample_subsets(factors::Array{OrderedDict{Symbol, Any}, 1},
                        ranges::Array{UnitRange{Int}, 1},
                        designs::Int;
                        check_bounds::Bool = true,
                        scale::Function = scale_boxdraper_encoding!,
                        compute_all_metrics::Bool = false)
    sampled_subsets = []

    for subset = 1:length(ranges)
        sampled_subset = sample_subset(factors[subset],
                                       ranges[subset],
                                       designs,
                                       check_bounds = check_bounds,
                                       scale = scale)

        push!(sampled_subsets, sampled_subset)
    end

    return sampled_subsets
end
