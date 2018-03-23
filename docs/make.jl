using Documenter, ExperimentalDesign

makedocs(format   = :html,
         sitename = "ExperimentalDesign",
         pages    = ["index.md"])

deploydocs(repo   = "github.com/phrb/ExperimentalDesign.jl.git",
           deps   = Deps.pip("mkdocs", "python-markdown-math"),
           julia  = "0.6")
