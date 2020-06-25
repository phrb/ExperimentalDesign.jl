"""
$(TYPEDEF)
"""
abstract type AbstractDesign end

"""
$(TYPEDEF)
"""
abstract type AbstractScreeningDesign <: AbstractDesign end

"""
$(TYPEDEF)
"""
abstract type AbstractFactorialDesign <: AbstractDesign end

"""
$(TYPEDEF)
"""
abstract type AbstractOptimalDesign <: AbstractDesign end

"""
$(TYPEDEF)
"""
abstract type AbstractRandomDesign <: AbstractDesign end

"""
$(TYPEDEF)

Encapsulates  a  Plackett-Burman  screening  design  constructed  using  Paley's
method.  Factor  levels are  encoded as  `:high` and  `:low` symbols,  and extra
dummy variables, possibly  aliasing interactions between factors,  will be added
to pad a design.

$(TYPEDFIELDS)
"""
struct PlackettBurman <: AbstractScreeningDesign
    matrix::DataFrame
    factors::Tuple
    dummy_factors::Tuple
    formula::FormulaTerm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> PlackettBurman(@formula(y ~ x1 + x2 + x3 + x4))
PlackettBurman
Dimension: (8, 7)
Factors: (:x1, :x2, :x3, :x4)
Dummy Factors: (:dummy1, :dummy2, :dummy3)
Formula: y ~ -1 + x1 + x2 + x3 + x4 + dummy1 + dummy2 + dummy3
Design Matrix:
8×7 DataFrame
│ Row │ x1    │ x2    │ x3    │ x4    │ dummy1 │ dummy2 │ dummy3 │
│     │ Int64 │ Int64 │ Int64 │ Int64 │ Int64  │ Int64  │ Int64  │
├─────┼───────┼───────┼───────┼───────┼────────┼────────┼────────┤
│ 1   │ 1     │ 1     │ 1     │ 1     │ 1      │ 1      │ 1      │
│ 2   │ -1    │ 1     │ -1    │ 1     │ 1      │ -1     │ -1     │
│ 3   │ 1     │ -1    │ 1     │ 1     │ -1     │ -1     │ -1     │
│ 4   │ -1    │ 1     │ 1     │ -1    │ -1     │ -1     │ 1      │
│ 5   │ 1     │ 1     │ -1    │ -1    │ -1     │ 1      │ -1     │
│ 6   │ 1     │ -1    │ -1    │ -1    │ 1      │ -1     │ 1      │
│ 7   │ -1    │ -1    │ -1    │ 1     │ -1     │ 1      │ 1      │
│ 8   │ -1    │ -1    │ 1     │ -1    │ 1      │ 1      │ -1     │

```
"""
function PlackettBurman(formula::FormulaTerm; symbol_encoding::Bool = false)
    symbol_factors = Tuple(r.sym for r in formula.rhs)

    initial_design = plackettburman(length(symbol_factors))
    categorical_design = similar(initial_design, Symbol)

    if symbol_encoding
        for i = 1:size(initial_design, 1), j = 1:size(initial_design, 2)
            categorical_design[i, j] = initial_design[i, j] == 1.0 ? :high : :low
        end

        design = DataFrame(categorical_design)
    else
        design = DataFrame(initial_design)
    end

    dummy_factors = Tuple(Symbol("dummy" * string(i)) for i = 1:(length(names(design)) -
                                                                 length(symbol_factors)))

    design_names = (symbol_factors..., dummy_factors...)
    rename!(design, collect(design_names))

    PlackettBurman(design,
                   symbol_factors,
                   dummy_factors,
                   term(formula.lhs) ~ term(-1) + sum(term.(design_names)))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> PlackettBurman(4)
PlackettBurman
Dimension: (8, 7)
Factors: (:factor1, :factor2, :factor3, :factor4)
Dummy Factors: (:dummy1, :dummy2, :dummy3)
Formula: 0 ~ -1 + factor1 + factor2 + factor3 + factor4 + dummy1 + dummy2 + dummy3
Design Matrix:
8×7 DataFrame
│ Row │ factor1 │ factor2 │ factor3 │ factor4 │ dummy1 │ dummy2 │ dummy3 │
│     │ Int64   │ Int64   │ Int64   │ Int64   │ Int64  │ Int64  │ Int64  │
├─────┼─────────┼─────────┼─────────┼─────────┼────────┼────────┼────────┤
│ 1   │ 1       │ 1       │ 1       │ 1       │ 1      │ 1      │ 1      │
│ 2   │ -1      │ 1       │ -1      │ 1       │ 1      │ -1     │ -1     │
│ 3   │ 1       │ -1      │ 1       │ 1       │ -1     │ -1     │ -1     │
│ 4   │ -1      │ 1       │ 1       │ -1      │ -1     │ -1     │ 1      │
│ 5   │ 1       │ 1       │ -1      │ -1      │ -1     │ 1      │ -1     │
│ 6   │ 1       │ -1      │ -1      │ -1      │ 1      │ -1     │ 1      │
│ 7   │ -1      │ -1      │ -1      │ 1       │ -1     │ 1      │ 1      │
│ 8   │ -1      │ -1      │ 1       │ -1      │ 1      │ 1      │ -1     │

```
"""
function PlackettBurman(factors::Int)
    PlackettBurman(ConstantTerm(0) ~ sum(term.(Symbol("factor" * string(i)) for i = 1:factors)))
end

"""
$(TYPEDEF)

Encapsulates  an *explicit*  full factorial  design, where  all experiments  are
completely computed  and stored.  Use [`IterableFullFactorial`](@ref)  to obtain
an *implicit* factorial  design with arbitrary access to the  elements of a full
factorial design.

$(TYPEDFIELDS)
"""
struct FullFactorial <: AbstractFactorialDesign
    matrix::DataFrame
    factors::NamedTuple
    formula::FormulaTerm
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct IterableFullFactorial <: AbstractFactorialDesign
    matrix::DataFrame
    factors::NamedTuple
    formula::FormulaTerm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> FullFactorial((A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0]), @formula(y ~ A + B +C))
FullFactorial
Dimension: (12, 3)
Factors: (A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0])
Formula: y ~ A + B + C
Design Matrix:
12×3 DataFrame
│ Row │ A   │ B   │ C    │
│     │ Any │ Any │ Any  │
├─────┼─────┼─────┼──────┤
│ 1   │ 1   │ a   │ 1.0  │
│ 2   │ 2   │ a   │ 1.0  │
│ 3   │ 4   │ a   │ 1.0  │
│ 4   │ 1   │ b   │ 1.0  │
│ 5   │ 2   │ b   │ 1.0  │
│ 6   │ 4   │ b   │ 1.0  │
│ 7   │ 1   │ a   │ -1.0 │
│ 8   │ 2   │ a   │ -1.0 │
│ 9   │ 4   │ a   │ -1.0 │
│ 10  │ 1   │ b   │ -1.0 │
│ 11  │ 2   │ b   │ -1.0 │
│ 12  │ 4   │ b   │ -1.0 │

```
"""
function FullFactorial(factors::NamedTuple, formula::FormulaTerm)
    iterator = fullfactorial(values(factors))
    matrix = DataFrame(explicit_fullfactorial(iterator))
    rename!(matrix, [r.sym for r in formula.rhs])
    FullFactorial(matrix, factors, formula)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> FullFactorial((A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0]))
FullFactorial
Dimension: (12, 3)
Factors: (A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0])
Formula: 0 ~ A + B + C
Design Matrix:
12×3 DataFrame
│ Row │ A   │ B   │ C    │
│     │ Any │ Any │ Any  │
├─────┼─────┼─────┼──────┤
│ 1   │ 1   │ a   │ 1.0  │
│ 2   │ 2   │ a   │ 1.0  │
│ 3   │ 4   │ a   │ 1.0  │
│ 4   │ 1   │ b   │ 1.0  │
│ 5   │ 2   │ b   │ 1.0  │
│ 6   │ 4   │ b   │ 1.0  │
│ 7   │ 1   │ a   │ -1.0 │
│ 8   │ 2   │ a   │ -1.0 │
│ 9   │ 4   │ a   │ -1.0 │
│ 10  │ 1   │ b   │ -1.0 │
│ 11  │ 2   │ b   │ -1.0 │
│ 12  │ 4   │ b   │ -1.0 │

```
"""
function FullFactorial(factors::NamedTuple)
    FullFactorial(factors, ConstantTerm(0) ~ sum(term.(keys(factors))))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> FullFactorial(([1, 2, 4], [:a, :b], [1.0, -1.0]))
FullFactorial
Dimension: (12, 3)
Factors: (factor1 = [1, 2, 4], factor2 = [:a, :b], factor3 = [1.0, -1.0])
Formula: 0 ~ factor1 + factor2 + factor3
Design Matrix:
12×3 DataFrame
│ Row │ factor1 │ factor2 │ factor3 │
│     │ Any     │ Any     │ Any     │
├─────┼─────────┼─────────┼─────────┤
│ 1   │ 1       │ a       │ 1.0     │
│ 2   │ 2       │ a       │ 1.0     │
│ 3   │ 4       │ a       │ 1.0     │
│ 4   │ 1       │ b       │ 1.0     │
│ 5   │ 2       │ b       │ 1.0     │
│ 6   │ 4       │ b       │ 1.0     │
│ 7   │ 1       │ a       │ -1.0    │
│ 8   │ 2       │ a       │ -1.0    │
│ 9   │ 4       │ a       │ -1.0    │
│ 10  │ 1       │ b       │ -1.0    │
│ 11  │ 2       │ b       │ -1.0    │
│ 12  │ 4       │ b       │ -1.0    │

```
"""
function FullFactorial(factors::Tuple)
    FullFactorial(NamedTuple{Tuple(Symbol("factor" * string(i)) for i = 1:length(factors))}(factors))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> FullFactorial(fill([-1, 1], 3))
FullFactorial
Dimension: (8, 3)
Factors: (factor1 = [-1, 1], factor2 = [-1, 1], factor3 = [-1, 1])
Formula: 0 ~ factor1 + factor2 + factor3
Design Matrix:
8×3 DataFrame
│ Row │ factor1 │ factor2 │ factor3 │
│     │ Int64   │ Int64   │ Int64   │
├─────┼─────────┼─────────┼─────────┤
│ 1   │ -1      │ -1      │ -1      │
│ 2   │ 1       │ -1      │ -1      │
│ 3   │ -1      │ 1       │ -1      │
│ 4   │ 1       │ 1       │ -1      │
│ 5   │ -1      │ -1      │ 1       │
│ 6   │ 1       │ -1      │ 1       │
│ 7   │ -1      │ 1       │ 1       │
│ 8   │ 1       │ 1       │ 1       │

```
"""
function FullFactorial(factors::Array)
    FullFactorial(tuple(factors...))
end


"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct FractionalFactorial <: AbstractFactorialDesign
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct RandomDesign <: AbstractRandomDesign
    matrix::DataFrame
    factors::NamedTuple
    formula::FormulaTerm
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

Encapsulates one or more `Distribution`s, generating random designs.  Receives a
`NamedTuple`  of  factor names  associated  with  `Distribution`s from  the  the
*Distributions* package.   After instantiating a `DesignDistribution`,  you must
request   samples   from   it    using   [`rand`](@ref),   which   generates   a
[`RandomDesign`](@ref)

"""
struct DesignDistribution
    factors::NamedTuple
    formula::FormulaTerm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> DesignDistribution((f1 = Uniform(2, 3), f2 = DiscreteUniform(-1, 5), f3 = Uniform(5, 10)))
DesignDistribution
Formula: 0 ~ f1 + f2 + f3
Factor Distributions:
f1: Uniform{Float64}(a=2.0, b=3.0)
f2: DiscreteUniform(a=-1, b=5)
f3: Uniform{Float64}(a=5.0, b=10.0)

```
"""
function DesignDistribution(factors::NamedTuple)
    DesignDistribution(factors, ConstantTerm(0) ~ sum(term.(keys(factors))))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> DesignDistribution((Uniform(2, 3), DiscreteUniform(-1, 5), Uniform(5, 10)))
DesignDistribution
Formula: 0 ~ factor1 + factor2 + factor3
Factor Distributions:
factor1: Uniform{Float64}(a=2.0, b=3.0)
factor2: DiscreteUniform(a=-1, b=5)
factor3: Uniform{Float64}(a=5.0, b=10.0)

```
"""
function DesignDistribution(factors::Tuple)
    DesignDistribution(NamedTuple{Tuple(Symbol("factor" * string(i)) for i = 1:length(factors))}(factors))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> DesignDistribution([Uniform(2, 3), DiscreteUniform(-1, 5), Uniform(5, 10)])
DesignDistribution
Formula: 0 ~ factor1 + factor2 + factor3
Factor Distributions:
factor1: Uniform{Float64}(a=2.0, b=3.0)
factor2: DiscreteUniform(a=-1, b=5)
factor3: Uniform{Float64}(a=5.0, b=10.0)

```
"""
function DesignDistribution(factors::Array)
    DesignDistribution(tuple(factors...))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> DesignDistribution(DiscreteNonParametric([-1, 1], [0.5, 0.5]), 6)
DesignDistribution
Formula: 0 ~ factor1 + factor2 + factor3 + factor4 + factor5 + factor6
Factor Distributions:
factor1: DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])
factor2: DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])
factor3: DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])
factor4: DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])
factor5: DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])
factor6: DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])

```
"""
function DesignDistribution(distribution::D, n::Int) where D <: Distribution
    factors = tuple(fill(distribution, n)...)
    DesignDistribution(factors)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> rand(DesignDistribution((f1 = Uniform(2, 3), f2 = DiscreteUniform(-1, 5), f3 = Uniform(5, 10))), 12)
ExperimentalDesign.RandomDesign
Dimension: (12, 3)
Factors: (f1 = Uniform{Float64}(a=2.0, b=3.0), f2 = DiscreteUniform(a=-1, b=5), f3 = Uniform{Float64}(a=5.0, b=10.0))
Formula: 0 ~ f1 + f2 + f3
Design Matrix:
12×3 DataFrame
│ Row │ f1      │ f2      │ f3      │
│     │ Float64 │ Float64 │ Float64 │
├─────┼─────────┼─────────┼─────────┤
│ 1   │ 2.04922 │ 5.0     │ 5.72029 │
│ 2   │ 2.59117 │ 1.0     │ 8.48192 │
│ 3   │ 2.77148 │ 4.0     │ 5.42517 │
│ 4   │ 2.25659 │ 5.0     │ 6.75815 │
│ 5   │ 2.64968 │ 5.0     │ 8.94993 │
│ 6   │ 2.31523 │ 4.0     │ 5.89413 │
│ 7   │ 2.86526 │ 3.0     │ 9.8807  │
│ 8   │ 2.44753 │ -1.0    │ 9.24986 │
│ 9   │ 2.08284 │ 3.0     │ 9.99911 │
│ 10  │ 2.81201 │ 2.0     │ 7.28848 │
│ 11  │ 2.64232 │ 4.0     │ 5.34507 │
│ 12  │ 2.08189 │ 4.0     │ 8.303   │

```
"""
function rand(distribution::DesignDistribution, n::Int = 1)
    r_design = zeros(n, length(distribution.factors))
    RandomDesign(DataFrame(random_design!(r_design,
                                          values(distribution.factors), n),
                           collect(keys(distribution.factors))),
                 distribution.factors,
                 distribution.formula)
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

Contains a  set of candidate experiments,  a target optimality criterion,  and a
model prior  formula.  A set  of experiments  that maximizes a  given optimality
criterion can selected from a candidate set using the [`kl_exchange`](@ref)
algorithm.

"""
struct OptimalDesign <: AbstractOptimalDesign
    matrix::DataFrame
    factors::NamedTuple
    formula::FormulaTerm
    selected_experiments::Array{Int}
    criteria::Dict{Symbol, Float64}
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> design_distribution = DesignDistribution((f1 = Uniform(2, 3), f2 = DiscreteUniform(-1, 5), f3 = Uniform(5, 10)))
DesignDistribution
Formula: 0 ~ f1 + f2 + f3
Factor Distributions:
f1: Uniform{Float64}(a=2.0, b=3.0)
f2: DiscreteUniform(a=-1, b=5)
f3: Uniform{Float64}(a=5.0, b=10.0)

julia> design = rand(design_distribution, 400);

julia> f = @formula 0 ~ f1 + f2 + f3 + f2 ^ 2;

julia> OptimalDesign(design, f, 10)
OptimalDesign
Dimension: (10, 3)
Factors: (f1 = Uniform{Float64}(a=2.0, b=3.0), f2 = DiscreteUniform(a=-1, b=5), f3 = Uniform{Float64}(a=5.0, b=10.0))
Formula: 0 ~ f1 + f2 + f3 + :(f2 ^ 2)
Selected Candidate Rows: [222, 103, 384, 140, 139, 156, 63, 184, 169, 54]
Optimality Criteria: Dict(:D => 2.418045723405036)
Design Matrix:
10×3 DataFrame
│ Row │ f1      │ f2      │ f3      │
│     │ Float64 │ Float64 │ Float64 │
├─────┼─────────┼─────────┼─────────┤
│ 1   │ 2.02407 │ -1.0    │ 5.07959 │
│ 2   │ 2.0203  │ 5.0     │ 9.50193 │
│ 3   │ 2.80676 │ -1.0    │ 9.87726 │
│ 4   │ 2.9996  │ 2.0     │ 5.03174 │
│ 5   │ 2.14846 │ 2.0     │ 9.72596 │
│ 6   │ 2.84392 │ 5.0     │ 5.01213 │
│ 7   │ 2.04046 │ 2.0     │ 5.00331 │
│ 8   │ 2.93432 │ 1.0     │ 9.71405 │
│ 9   │ 2.79858 │ -1.0    │ 5.35685 │
│ 10  │ 2.97635 │ 5.0     │ 8.10695 │

```
"""
function OptimalDesign(candidates::AbstractDesign,
                       formula::FormulaTerm,
                       experiments::Int;
                       tolerance::Float64 = 1e-9,
                       seed_design_size::Int = 2,
                       max_iterations::Int = 1000,
                       design_k::Int = experiments,
                       candidates_l::Int = size(candidates.matrix, 1) - design_k)
    selected = kl_exchange(formula,
                           candidates.matrix,
                           experiments = experiments,
                           tolerance = tolerance,
                           seed_design_size = seed_design_size,
                           max_iterations = max_iterations,
                           design_k = design_k,
                           candidates_l = candidates_l)

    optimality = d_criterion(selected)

    OptimalDesign(candidates.matrix[selected.indices[1], :],
                  candidates.factors,
                  formula,
                  selected.indices[1],
                  Dict(:D => optimality))
end
