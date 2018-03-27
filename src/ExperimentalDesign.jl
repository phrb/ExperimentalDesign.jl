#__precompile__()

module ExperimentalDesign

using Reexport

@reexport using DataFrames
@reexport using DataStructures

using IterTools,
      StatsBase,
      StatPlots,
      Primes

export plackett_burman,
       is_plackett_burman,
       paley

export plot_subsets,
       sample_subsets,
       sample_subset,
       generate_designs

export scale_orthogonal!,
       scale_boxdraper_encoding!,
       generate_model_matrix,
       get_prediction_variances,
       build_linear_formula,
       get_model_variables

export condition_number,
       d_optimality,
       a_optimality,
       v_optimality,
       g_optimality,
       g_efficiency,
       d_efficiency_lower_bound,
       d_efficiency_lower_bound_algdesign

include("plackett_burman.jl")
include("d_optimal.jl")

end
