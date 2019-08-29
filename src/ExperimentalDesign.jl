#__precompile__()

module ExperimentalDesign

using Random
using Logging
using Primes
using DataFrames
using Distributions
using StatsModels
using DocStringExtensions
using LinearAlgebra

import Random: rand

# Types

export    AbstractDesign,   AbstractScreeningDesign,    AbstractFactorialDesign,
    AbstractOptimalDesign,  PlackettBurman, FullFactorial,  FractionalFactorial,
    OptimalDesign, RandomDesign

# Methods

export fullfactorial,  explicit_fullfactorial, plackettburman, isplackettburman,
    paley,     optimize_design,     expanded_design,     generate_random_design,
    rand, random_design, d_criterion, next_offset_divisible_prime

include("design.jl")
include("plackett_burman.jl")
include("factorial.jl")
include("random.jl")
include("d_optimal/variance_predictions.jl")

end
