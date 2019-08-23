#__precompile__()

module ExperimentalDesign

using Random
using Logging
using Primes
using DataFrames
using StatsModels
using DocStringExtensions
using LinearAlgebra

# Plackett-Burman Designs

export plackettburman, isplackettburman, paley

# D-Optimal Designs

export optimize_design, expanded_design, random_design, d_criterion

include("design.jl")
include("plackett_burman.jl")
include("d_optimal/variance_predictions.jl")

end
