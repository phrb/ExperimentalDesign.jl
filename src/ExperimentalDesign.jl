#__precompile__()

module ExperimentalDesign

using Random
using Combinatorics
using Logging
using Primes
using DataFrames
using Distributions
using StatsModels
using DocStringExtensions
using LinearAlgebra
using LatinHypercubeSampling

import Random: rand

# Types

export    AbstractDesign,   AbstractScreeningDesign,    AbstractFactorialDesign,
    AbstractOptimalDesign,  BoxBehnken, CentralComposite, PlackettBurman, FullFactorial, FractionalFactorial,
    FractionalFactorial2Level, OptimalDesign, DesignDistribution, CategoricalFactor,
    OptimLHCDesign, RandomLHCDesign

# Methods

export fullfactorial,  explicit_fullfactorial, plackettburman, fold!, isplackettburman,
    paley, rand, random_design!, next_offset_divisible_prime, kl_exchange, d_criterion,
    boxbehnken, ccdesign

# Pretty printing

export show

include("design.jl")
include("plackett_burman.jl")
include("factorial.jl")
include("random.jl")
include("kl_exchange.jl")
include("custom_show.jl")
include("categorical.jl")
include("central_composite.jl")
include("fold.jl")
include("box-behnken.jl")

end
