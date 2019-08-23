using Test
using ExperimentalDesign

tests = ["variance_predictions.jl"]

@testset "ExperimentalDesign" begin
    for test in tests
        include(test)
    end
end
