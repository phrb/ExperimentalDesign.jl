using Test
using Documenter
using ExperimentalDesign

tests = []

@testset "ExperimentalDesign" begin
    for test in tests
        include(test)
    end

    @testset "Doctests" begin
        DocMeta.setdocmeta!(ExperimentalDesign,
                            :DocTestSetup,
                            :(using ExperimentalDesign, Distributions,
                              Random, StatsModels, DataFrames;
                              Random.seed!(443591););
                            recursive = true)
        doctest(ExperimentalDesign)
    end
end
