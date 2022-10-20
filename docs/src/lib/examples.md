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

p = @df design.matrix scatter(:size,
    :weight,
    size = (600, 600),
    xlabel = "size",
    ylabel = "weight",
    xlim = [23.0, 32.0],
    ylim = [0.0, 100.0],
    legend = false,
    title = "Uniformly Sampled Design")

png(p, "plot1.png")
nothing
```

![](plot1.png)

#### Generating Experiments for a Linear Hypothesis

```@example edexample
design = rand(design_distribution, 400);

f = @formula 0 ~ size + weight

optimal_design = OptimalDesign(design, f, 10)

p = @df optimal_design.matrix scatter(:size,
    :weight,
    size = (600, 600),
    xlabel = "size",
    ylabel = "weight",
    xlim = [23.0, 32.0],
    ylim = [0.0, 100.0],
    legend = false,
    title = "Design for y = size + weight")
png(p, "plot2.png")
nothing
```

![](plot2.png)

#### Generating Experiments for Other Terms


```@example edexample


design = rand(design_distribution, 400);
f = @formula 0 ~ size + weight + size ^ 2 + (1 / weight)

optimal_design = OptimalDesign(design, f, 20)

p = @df optimal_design.matrix scatter(:size,
    :weight,
    size = (600, 600),
    xlabel = "size",
    ylabel = "weight",
    xlim = [23.0, 32.0],
    ylim = [0.0, 100.0],
    legend = false,
    title = "Design for y = size + weight + (size ^ 2) + (1 / weight)")
png(p, "plot3.png")
nothing
```

![](plot3.png)

```@example edexample

design = rand(design_distribution, 800);
f = @formula 0 ~ size + weight + size ^ 2
optimal_design = OptimalDesign(design, f, 10)

p = @df optimal_design.matrix scatter(:size,
    :weight,
    size = (600, 600),
    xlabel = "size",
    ylabel = "weight",
    legend = false,
    title = "Design for y = size + weight + (size ^ 2)")

png(p, "plot4.png")
nothing
```

![](plot4.png)

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
p = @df optimal_design.matrix scatter(:f1,
    :f2,
    size = (600, 600),
    xlabel = "f1",
    ylabel = "f2",
    legend = false,
    title = "Optimal Design for y = f1 + f2")
png(p, "plot5.png")
nothing
```

![](plot5.png)


### Screening with Plackett-Burman Designs

#### Generating Plackett-Burman Designs

A Plackett-Burman design is an orthogonal design matrix for factors $f_1,\dots,f_N$. Factors are encoded by high and low values, which can be mapped to the interval $[-1, 1]$. For designs in this package, the design matrix is a `DataFrame` from the [DataFrame package](https://juliastats.org/GLM.jl/stable/). For example, let's create a Plackett-Burman design for 6 factors:

```@example edexample
design = PlackettBurman(6)
design.matrix
```

Note that it is not possible to construct exact Plackett-Burman designs for all numbers of factors. In the example above, we needed a seventh extra "dummy" column to construct the design for six factors.

Using the `PlackettBurman` constructor enables quick construction of minimal screening designs for scenarios where we ignore interactions. We can access the underlying formula, which is a `Term` object from the [StatsModels package](https://juliastats.org/StatsModels.jl/stable/):


```@example edexample
println(design.formula)
```

Notice we ignore interactions and include the dummy factor in the model. Strong main effects attributed to dummy factors may indicate important interactions.

We can obtain a tuple with the names of dummy factors:


```@example edexample
design.dummy_factors
```

We can also get the main factors tuple:

```@example edexample
design.factors
```

You can check other constructors on [the docs](https://phrb.github.io/ExperimentalDesign.jl/dev/lib/public/#ExperimentalDesign.PlackettBurman-Tuple{Int64}).

#### Computing Main Effects

Suppose that the response variable on the experiments specified in our screening design is computed by:

$$
y = 1.2 + (2.3f_1) + (-3.4f_2) + (7.12f_3) + (-0.03f_4) + (1.1f_5) + (-0.5f_6) + \varepsilon
$$

The coefficients we want to estimate are:

| Intercept | factor1 | factor2 | factor3 | factor4 | factor5 | factor6 |
|---|---|---|---|---|---|---|
| 1.2 | 2.3 | -3.4 | 7.12 | -0.03 | 1.1 | -0.5 |

The corresponding Julia function is:

```@example edexample
function y(x)
    return (1.2) +
           (2.3 * x[1]) +
           (-3.4 * x[2]) +
           (7.12 * x[3]) +
           (-0.03 * x[4]) +
           (1.1 * x[5]) +
           (-0.5 * x[6]) +
           (1.1 * randn())
end
```

We can compute the response column for our design using the cell below. Recall that the default is to call the response column `:response`. We are going to set the seeds each time we run `y(x)`, so we analyse same results. Play with different seeds to observe variability of estimates.

```@example edexample
Random.seed!(192938)

design.matrix[!, :response] = y.(eachrow(design.matrix[:, collect(design.factors)]))
design.matrix
```

Now, we use the `lm` function from the [GLM package](https://juliastats.org/GLM.jl/stable/) to fit a linear model using the design's matrix and formula:

```@example edexample
lm(term(:response) ~ design.formula.rhs, design.matrix)
```

The table below shows the coefficients estimated by the linear model fit using the Plackett-Burman Design. The purpose of a screening design is not to estimate the actual coefficients, but instead to compute factor main effects. Note that standard errors are the same for every factor estimate. This happens because the design is orthogonal.

| | Intercept | factor1 | factor2 | factor3 | factor4 | factor5 | factor6 | dummy1 |
|---|---|---|---|---|---|---|---|---|
| Original | 1.2 | 2.3 | -3.4 | 7.12 | -0.03 | 1.1 | -0.5 | $-$ |
| Plackett-Burman Main Effects | $-$ | 2.66359 | -2.80249 | 7.71644 | -0.0497774 | 0.612681 | -0.112675 | -0.437176 |

We can use the coefficient magnitudes to infer that factor 3 probably has a strong main effect, and that factor 6 has not. Our dummy column had a relatively small coefficient estimate, so we could attempt to ignore interactions on subsequent experiments.

#### Fitting a Linear Model

We can also try to fit a linear model on our design data in order to estimate coefficients. We would need to drop the dummy column and add the intercept term:

```@example edexample
lm(term(:response) ~ sum(term.(design.factors)), design.matrix)
```

Our table so far looks like this:

| | Intercept | factor1 | factor2 | factor3 | factor4 | factor5 | factor6 | dummy1 |
|---|---|---|---|---|---|---|---|---|
| Original | 1.2 | 2.3 | -3.4 | 7.12 | -0.03 | 1.1 | -0.5 | $-$ |
| Plackett-Burman Main Effects | $-$ | 2.66359 | -2.80249 | 7.71644 | -0.0497774 | 0.612681 | -0.112675 | -0.437176 |
| Plackett-Burman Estimate | 0.940183 | 2.66359 | -2.80249 | 7.71644 | -0.0497774 | 0.612681 | -0.112675 | $-$ |

Notice that, since the standard errors are the same for all factors, factors with stronger main effects are better estimated. Notice that, despite the "good" coefficient estimates, the confidence intervals are really large.

This is a biased comparison where the screening design "works" for coefficient estimation as well, but we would rather use fractional factorial or optimal designs to estimate the coefficients of factors with strong effects. Screening should be used to compute main effects and identifying which factors to test next.

#### Generating Random Designs

We can also compare the coefficients produced by the same linear model fit, but using a random design. For more information, check [the docs](https://phrb.github.io/ExperimentalDesign.jl/dev/lib/public/#ExperimentalDesign.RandomDesign-Tuple{NamedTuple}).

```@example edexample
Random.seed!(8418172)

design_distribution = DesignDistribution(DiscreteNonParametric([-1, 1], [0.5, 0.5]), 6)
random_design = rand(design_distribution, 8)

random_design.matrix[!, :response] = y.(eachrow(random_design.matrix[:, :]))
random_design.matrix
```

```@example edexample
lm(term(:response) ~ random_design.formula.rhs, random_design.matrix)
```

Now, our table looks like this:

| | Intercept | factor1 | factor2 | factor3 | factor4 | factor5 | factor6 | dummy1 |
|---|---|---|---|---|---|---|---|---|
| Original | 1.2 | 2.3 | -3.4 | 7.12 | -0.03 | 1.1 | -0.5 | $-$ |
| Plackett-Burman Main Effects | $-$ | 2.66359 | -2.80249 | 7.71644 | -0.0497774 | 0.612681 | -0.112675 | -0.437176 |
| Plackett-Burman Estimate | 0.940183 | 2.66359 | -2.80249 | 7.71644 | -0.0497774 | 0.612681 | -0.112675 | $-$ |
| Single Random Design Estimate | 0.761531 | 1.67467 | -3.05027 | 7.72484 | 0.141204 | 1.71071 | -0.558869 | $-$ |

The estimates produced using random designs will have larger confidence intervals, and therefore increased variability. The Plackett-Burman design is fixed, but can be randomised. The variability of main effects estimates using screening designs will depend on measurement or model error.

#### Generating Full Factorial Designs

In this toy example, it is possible to generate all the possible combinations of six binary factors and compute the response. Although it costs 64 experiments, the linear model fit for the full factorial design should produce the best coefficient estimates.

The simplest full factorial design constructor receives an array of possible factor levels. For more, check [the docs](https://phrb.github.io/ExperimentalDesign.jl/dev/lib/public/#ExperimentalDesign.FullFactorial-Tuple{NamedTuple,StatsModels.FormulaTerm}).

```@example edexample
Random.seed!(2989476)

factorial_design = FullFactorial(fill([-1, 1], 6))
factorial_design.matrix[!, :response] = y.(eachrow(factorial_design.matrix[:, :]))

lm(term(:response) ~ factorial_design.formula.rhs, factorial_design.matrix)
```

The confidence intervals for this fit are much smaller. Since we have all information on all factors and this is a balanced design, the standard error is the same for all estimates. Here's the complete table:

| | Intercept | factor1 | factor2 | factor3 | factor4 | factor5 | factor6 | dummy1 |
|---|---|---|---|---|---|---|---|---|
| Original | 1.2 | 2.3 | -3.4 | 7.12 | -0.03 | 1.1 | -0.5 | $-$ |
| Plackett-Burman Main Effects | $-$ | 2.66359 | -2.80249 | 7.71644 | -0.0497774 | 0.612681 | -0.112675 | -0.437176 |
| Plackett-Burman Estimate | 0.940183 | 2.66359 | -2.80249 | 7.71644 | -0.0497774 | 0.612681 | -0.112675 | $-$ |
| Single Random Design Estimate | 0.600392 | 2.2371 | -2.56857 | 8.05743 | 0.140622 | 0.907918 | -0.600354 | $-$ |
| Full Factorial Estimate | 1.13095 | 2.23668 | -3.4775 | 6.95531 | -0.160546 | 0.975471 | -0.357748 | $-$ |

Full factorial designs may be too expensive in actual applications. Fractional factorial designs or optimal designs can be used to decrease costs while still providing good estimates. Screening designs are extremely cheap, and can help determine which factors can potentially be dropped on more expensive and precise designs.
