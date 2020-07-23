using Test,
    Documenter,
    ExperimentalDesign,
    Random,
    Distributions,
    StatsModels

# tests = ["kl_exchange.jl", "doctests.jl"]
tests = ["kl_exchange.jl"]

@testset "ExperimentalDesign" begin
    for test in tests
        include(test)
    end
end
