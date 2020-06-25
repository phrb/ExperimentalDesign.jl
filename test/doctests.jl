@testset "Doctests" begin
    random_seed = 443591
    DocMeta.setdocmeta!(ExperimentalDesign,
                        :DocTestSetup,
                        :(using ExperimentalDesign, Distributions,
                          Random, StatsModels, DataFrames;
                          Random.seed!($random_seed););
                        recursive = true)
    doctest(ExperimentalDesign)
end
