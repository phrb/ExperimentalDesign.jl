function replace_zero(x, tol = 1e-4)
    return isapprox(x, 0.0, atol = tol) ? 0.0 : x
end

function replace_inf(x, lim = 1e5)
    if x == Inf
        return lim
    elseif x == -Inf
        return -lim
    else
        return x
    end
end

function plot_subsets(sampled_subsets, factor_title = " "; columns = [:D, :DELB])
    upscale    = 2
    small_font = Plots.font("sans-serif", 10.0 * upscale)
    large_font = Plots.font("sans-serif", 14.0 * upscale)
    default(titlefont  = large_font,
            guidefont  = large_font,
            tickfont   = small_font,
            legendfont = small_font)
    default(size = (700 * upscale, 900 * upscale))
    default(dpi = 300)

    plotly()

    subplots = []

    global_max_d = 0.0

    sampled_subsets = deepcopy(sampled_subsets)

    for subset in sampled_subsets
        global_max_d = max(global_max_d, subset[:D]...)
    end

    global_max_x = global_max_d == 0.0 ? string(global_max_d) : string("1e",
                                                                       floor(log(10, global_max_d)))

    for subset in sampled_subsets
        for column in columns
            subset[column] = replace_inf(subset[column])
            subset[column] = replace_zero.(subset[column])
        end

        max_d = max(subset[:D]...)
        max_x = max_d == 0.0 ? string(max_d) : string("1e",
                                                      floor(log(10, max_d)))

        push!(subplots,
              histogram(Array(subset[:D]), labels = "Designs",
                        normalize = :probability,
                        title  = string("D-Optimality for ", subset[:Length][1],
                                        " Experiments, ", size(subset, 1),
                                        " Samples and ", factor_title, " Factors"),
                        xticks = ([0, max_d, global_max_d],
                                  ["0", max_x, global_max_x]),
                        xlims = (0.0, global_max_d),
                        color  = :lightblue),
              histogram(Array(subset[:DELB]), labels = "Designs",
                        normalize = :probability,
                        title = string("D-Efficiency for ", subset[:Length][1],
                                        " Experiments, ", size(subset, 1),
                                        " Samples and ", factor_title, " Factors"),
                        xlims = (0.0, 1.0),
                        color = :lightgreen))
    end

    plot(subplots..., legend = false,
         layout = (length(sampled_subsets), 2))
end
