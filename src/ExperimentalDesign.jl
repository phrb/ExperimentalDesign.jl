#__precompile__()

module ExperimentalDesign

using Random, Logging, Primes, DataFrames, StatsModels, DocStringExtensions, LinearAlgebra

# Types

export    AbstractDesign,   AbstractScreeningDesign,    AbstractFactorialDesign,
    AbstractOptimalDesign,  PlackettBurman, FullFactorial,  FractionalFactorial,
    Optimal

# Methods

export fullfactorial,  explicit_fullfactorial, plackettburman, isplackettburman,
    paley,   optimize_design,   expanded_design,   random_design,   d_criterion,
    next_offset_divisible_prime

include("design.jl")
include("plackett_burman.jl")
include("factorial.jl")
include("d_optimal/variance_predictions.jl")

end
