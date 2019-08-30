# ExperimentalDesign.jl

| Build | Docs | Test Coverage |
| --- | --- | --- |
| [![Build Status](https://travis-ci.org/phrb/ExperimentalDesign.jl.svg?branch=master)](https://travis-ci.org/phrb/ExperimentalDesign.jl) | [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://phrb.github.io/ExperimentalDesign.jl/dev) | [![Coverage Status](https://coveralls.io/repos/github/phrb/ExperimentalDesign.jl/badge.svg?branch=master)](https://coveralls.io/github/phrb/ExperimentalDesign.jl?branch=master) [![codecov.io](http://codecov.io/github/phrb/ExperimentalDesign.jl/coverage.svg?branch=master)](http://codecov.io/github/phrb/ExperimentalDesign.jl?branch=master) |

ExperimentalDesign.jl  provides  tools  for  Design  of  Experiments  in  Julia,
enabling the construction  of designs for screening,  modeling, exploration, and
optimization.

Development  on this  package  is  ongoing, so  expect  things  to change.  Pull
requests are more than welcome!  Current features are:

- Designs supporting categorical and continuous factors
- Integration with [StatsModels](https://github.com/JuliaStats/StatsModels.jl) `@formula`
- Full factorial designs
- Plackett-Burman designs for screening
- Flexible random designs using the [Distributions](https://github.com/JuliaStats/Distributions.jl) package
- Variance-optimizing designs for several criteria

Intended features include the ones provided by R packages such as
[DoE.base](https://cran.r-project.org/web/packages/DoE.base/index.html),
[FrF2](https://cran.r-project.org/web/packages/FrF2/index.html), and
[AlgDesign](https://cran.r-project.org/web/packages/AlgDesign/index.html).

## Library Outline

```@contents
Pages = ["lib/public.md", "lib/internals.md"]
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md", "lib/internals.md"]
```
