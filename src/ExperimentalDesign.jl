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
    paley, rand, random_design, next_offset_divisible_prime

# Pretty printing

export show

include("design.jl")
include("plackett_burman.jl")
include("factorial.jl")
include("random.jl")
include("custom_show.jl")

end
