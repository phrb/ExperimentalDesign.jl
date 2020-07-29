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
    OptimalDesign, DesignDistribution, CategoricalFactor

# Methods

export fullfactorial,  explicit_fullfactorial, plackettburman, isplackettburman,
    paley, rand, random_design!, next_offset_divisible_prime, kl_exchange, d_criterion

# Pretty printing

export show

include("design.jl")
include("plackett_burman.jl")
include("factorial.jl")
include("random.jl")
include("kl_exchange.jl")
include("custom_show.jl")
include("categorical.jl")

end
