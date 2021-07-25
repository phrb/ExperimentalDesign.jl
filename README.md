# ExperimentalDesign

| Build | Docs | Test Coverage |
| --- | --- | --- |
| [![CI](https://github.com/phrb/ExperimentalDesign.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/phrb/ExperimentalDesign.jl/actions/workflows/CI.yml) | [![Coverage Status](https://coveralls.io/repos/github/phrb/ExperimentalDesign.jl/badge.svg?branch=master)](https://coveralls.io/github/phrb/ExperimentalDesign.jl?branch=master) [![codecov.io](http://codecov.io/github/phrb/ExperimentalDesign.jl/coverage.svg?branch=master)](http://codecov.io/github/phrb/ExperimentalDesign.jl?branch=master) |

ExperimentalDesign  provides  tools  for  Design  of  Experiments  in  Julia,
enabling the construction  of designs for screening,  modeling, exploration, and
optimization.

Development  on this  package  is  ongoing, so  expect  things  to change.  Pull
requests are more than welcome!

Check the [documentation](https://phrb.github.io/ExperimentalDesign.jl/dev/)
for the latest features and API, and check the examples directory for
Jupyter Notebooks and code.

Current features are:

- Designs that support categorical and continuous factors
- Integration with [StatsModels](https://github.com/JuliaStats/StatsModels.jl) `@formula`
- Full factorial designs:
  - Explicit: for small designs that fit in memory
  - Iterable: for larger designs, generates experiments on demand
- Plackett-Burman designs for screening (check the [example](https://github.com/phrb/ExperimentalDesign.jl/blob/master/examples/Screening%20with%20Plackett-Burman%20Designs.ipynb))
- Flexible random designs using the [Distributions](https://github.com/JuliaStats/Distributions.jl) package
- Several variance-optimizing criteria

Intended features include the ones provided by R packages such as
[DoE.base](https://cran.r-project.org/web/packages/DoE.base/index.html),
[FrF2](https://cran.r-project.org/web/packages/FrF2/index.html), and
[AlgDesign](https://cran.r-project.org/web/packages/AlgDesign/index.html).
