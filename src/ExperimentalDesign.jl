#__precompile__()

module ExperimentalDesign

using StatPlots

import StatsBase: counts

import Primes: primes

import DataStructures: OrderedDict

import DataFrames: DataFrame,
                   rename!

import StatsModels: Formula,
                    Terms

# Plackett-Burman Designs

export plackett_burman,
       is_plackett_burman,
       paley

# D-Optimal Designs

export plot_subsets,
       sample_subsets,
       sample_subset,
       full_factorial_subset,
       get_expanded_values,
       expand_design,
       expand_factors,
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

include("plackett_burman/plackett_burman.jl")
include("d_optimal/modeling.jl")
include("d_optimal/criteria.jl")
include("d_optimal/sampling.jl")
include("d_optimal/plotting.jl")

end
