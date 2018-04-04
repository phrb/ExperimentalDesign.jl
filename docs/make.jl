using Documenter, ExperimentalDesign

makedocs(format   = :html,
         sitename = "ExperimentalDesign",
         pages    = ["index.md"])

deploydocs(repo   = "github.com/phrb/ExperimentalDesign.jl.git",
           target = "build",
           make   = nothing,
           deps   = nothing,
           julia  = "0.6.2")
