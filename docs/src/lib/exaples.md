### Optimal Design with KL-Exchange

#### Generating Random Designs (Uniform)

```@example edexample
using ExperimentalDesign, StatsModels, GLM, DataFrames, Distributions, Random, StatsPlots

design_distribution = DesignDistribution((size = Uniform(23, 32), weight = Uniform(0, 100)))
```

```@example edexample

rand(design_distribution, 3)
```

```@example edexample

design = rand(design_distribution, 400)

@df design.matrix scatter(:size,
    :weight,
    size = (600, 600),
    xlabel = "size",
    ylabel = "weight",
    xlim = [23.0, 32.0],
    ylim = [0.0, 100.0],
    legend = false,
    title = "Uniformly Sampled Design")
```

#### Designs with Categorical Factors


```@example edexample
design_distribution = DesignDistribution((f1 = DiscreteUniform(0, 5),
        f2 = CategoricalFactor(["cf", "cg", "ca"])))

```

```@example edexample
design = rand(design_distribution, 300);
f = @formula 0 ~ f1 + f1 ^ 2 + f2

optimal_design = OptimalDesign(design, f, 10)
```

```@example edexample
@df optimal_design.matrix scatter(:f1,
    :f2,
    size = (600, 600),
    xlabel = "f1",
    ylabel = "f2",
    legend = false,
    title = "Optimal Design for y = f1 + f2")
```
