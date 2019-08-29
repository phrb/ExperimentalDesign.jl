using Documenter, ExperimentalDesign

# The DOCSARGS environment variable can be used to pass additional arguments to make.jl.
# This is useful on CI, if you need to change the behavior of the build slightly but you
# can not change the .travis.yml or make.jl scripts any more (e.g. for a tag build).
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        push!(ARGS, arg)
    end
end

DocMeta.setdocmeta!(ExperimentalDesign,
                    :DocTestSetup,
                    :(using ExperimentalDesign, Distributions,
                      Random, StatsModels, DataFrames;
                      Random.seed!(443591););
                    recursive = true)

makedocs(
    format = Documenter.HTML(prettyurls = !("local" in ARGS)),
    modules = [ExperimentalDesign],
    sitename = "ExperimentalDesign.jl",
    authors = "Pedro Bruel and contributors.",
    linkcheck = !("skiplinks" in ARGS),
    strict = !("strict=false" in ARGS),
    doctest = ("doctest=only" in ARGS) ? :only : true,
    pages = [
        "Home" => "index.md",
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md"
        ],
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(repo = "github.com/phrb/ExperimentalDesign.jl.git")
