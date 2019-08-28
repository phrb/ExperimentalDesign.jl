using Test, Documenter, ExperimentalDesign

tests = ["variance_predictions.jl"]

@testset "ExperimentalDesign" begin
    for test in tests
        include(test)
    end

    @testset "Doctests" begin
        DocMeta.setdocmeta!(ExperimentalDesign,
                            :DocTestSetup,
                            :(using ExperimentalDesign, StatsModels, DataFrames);
                            recursive = true)
        doctest(ExperimentalDesign)
    end
end
