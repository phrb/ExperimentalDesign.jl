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
 Row │ x1     x2     x3     x4     dummy1  dummy2  dummy3
     │ Int64  Int64  Int64  Int64  Int64   Int64   Int64
─────┼────────────────────────────────────────────────────
   1 │     1      1      1      1       1       1       1
   2 │    -1      1     -1      1       1      -1      -1
   3 │     1     -1      1      1      -1      -1      -1
   4 │    -1      1      1     -1      -1      -1       1
   5 │     1      1     -1     -1      -1       1      -1
   6 │     1     -1     -1     -1       1      -1       1
   7 │    -1     -1     -1      1      -1       1       1
   8 │    -1     -1      1     -1       1       1      -1

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

        design = DataFrame(categorical_design, :auto)
    else
        design = DataFrame(initial_design, :auto)
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
julia> fold!(PlackettBurman(2))
8×3 DataFrame
│ Row │ factor1 │ factor2 │ dummy1 │
│     │ Int64   │ Int64   │ Int64  │
├─────┼─────────┼─────────┼────────┤
│ 1   │ 1       │ 1       │ 1      │
│ 2   │ 1       │ -1      │ -1     │
│ 3   │ -1      │ -1      │ 1      │
│ 4   │ -1      │ 1       │ -1     │
│ 5   │ -1      │ -1      │ -1     │
│ 6   │ -1      │ 1       │ 1      │
│ 7   │ 1       │ 1       │ -1     │
│ 8   │ 1       │ -1      │ 1      │

```
"""
function fold!(design::AbstractScreeningDesign)
    mirror_matrix = deepcopy(design.matrix)
    if mirror_matrix[1,1] isa Symbol
        for col in names(mirror_matrix)
            replace!(mirror_matrix[col], :low => :high, :high => :low)
        end
    else
        for col in names(mirror_matrix)
            replace!(mirror_matrix[col], 1 => -1, -1 => 1)
        end
    end
    append!(design.matrix, mirror_matrix)
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
 Row │ factor1  factor2  factor3  factor4  dummy1  dummy2  dummy3
     │ Int64    Int64    Int64    Int64    Int64   Int64   Int64
─────┼────────────────────────────────────────────────────────────
   1 │       1        1        1        1       1       1       1
   2 │      -1        1       -1        1       1      -1      -1
   3 │       1       -1        1        1      -1      -1      -1
   4 │      -1        1        1       -1      -1      -1       1
   5 │       1        1       -1       -1      -1       1      -1
   6 │       1       -1       -1       -1       1      -1       1
   7 │      -1       -1       -1        1      -1       1       1
   8 │      -1       -1        1       -1       1       1      -1

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
 Row │ A    B    C
     │ Any  Any  Any
─────┼────────────────
   1 │ 1    a    1.0
   2 │ 2    a    1.0
   3 │ 4    a    1.0
   4 │ 1    b    1.0
   5 │ 2    b    1.0
   6 │ 4    b    1.0
   7 │ 1    a    -1.0
   8 │ 2    a    -1.0
   9 │ 4    a    -1.0
  10 │ 1    b    -1.0
  11 │ 2    b    -1.0
  12 │ 4    b    -1.0

```
"""
function FullFactorial(factors::NamedTuple, formula::FormulaTerm)
    iterator = fullfactorial(values(factors))
    matrix = DataFrame(explicit_fullfactorial(iterator), :auto)
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
 Row │ A    B    C
     │ Any  Any  Any
─────┼────────────────
   1 │ 1    a    1.0
   2 │ 2    a    1.0
   3 │ 4    a    1.0
   4 │ 1    b    1.0
   5 │ 2    b    1.0
   6 │ 4    b    1.0
   7 │ 1    a    -1.0
   8 │ 2    a    -1.0
   9 │ 4    a    -1.0
  10 │ 1    b    -1.0
  11 │ 2    b    -1.0
  12 │ 4    b    -1.0

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
 Row │ factor1  factor2  factor3
     │ Any      Any      Any
─────┼───────────────────────────
   1 │ 1        a        1.0
   2 │ 2        a        1.0
   3 │ 4        a        1.0
   4 │ 1        b        1.0
   5 │ 2        b        1.0
   6 │ 4        b        1.0
   7 │ 1        a        -1.0
   8 │ 2        a        -1.0
   9 │ 4        a        -1.0
  10 │ 1        b        -1.0
  11 │ 2        b        -1.0
  12 │ 4        b        -1.0

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
 Row │ factor1  factor2  factor3
     │ Int64    Int64    Int64
─────┼───────────────────────────
   1 │      -1       -1       -1
   2 │       1       -1       -1
   3 │      -1        1       -1
   4 │       1        1       -1
   5 │      -1       -1        1
   6 │       1       -1        1
   7 │      -1        1        1
   8 │       1        1        1

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
f1: Distributions.Uniform{Float64}(a=2.0, b=3.0)
f2: Distributions.DiscreteUniform(a=-1, b=5)
f3: Distributions.Uniform{Float64}(a=5.0, b=10.0)

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
factor1: Distributions.Uniform{Float64}(a=2.0, b=3.0)
factor2: Distributions.DiscreteUniform(a=-1, b=5)
factor3: Distributions.Uniform{Float64}(a=5.0, b=10.0)

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
factor1: Distributions.Uniform{Float64}(a=2.0, b=3.0)
factor2: Distributions.DiscreteUniform(a=-1, b=5)
factor3: Distributions.Uniform{Float64}(a=5.0, b=10.0)

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
factor1: Distributions.DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor2: Distributions.DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor3: Distributions.DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor4: Distributions.DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor5: Distributions.DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor6: Distributions.DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])

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
Factors: (f1 = Distributions.Uniform{Float64}(a=2.0, b=3.0), f2 = Distributions.DiscreteUniform(a=-1, b=5), f3 = Distributions.Uniform{Float64}(a=5.0, b=10.0))
Formula: 0 ~ f1 + f2 + f3
Design Matrix:
12×3 DataFrame
 Row │ f1       f2       f3
     │ Float64  Float64  Float64
─────┼───────────────────────────
   1 │ 2.04922     -1.0  9.62061
   2 │ 2.59117      3.0  5.79254
   3 │ 2.77148     -1.0  6.22902
   4 │ 2.25659      4.0  6.00256
   5 │ 2.64968      2.0  9.21818
   6 │ 2.31523      5.0  8.80749
   7 │ 2.86526      2.0  8.47297
   8 │ 2.44753      0.0  6.11722
   9 │ 2.08284     -1.0  6.34626
  10 │ 2.81201     -1.0  6.44508
  11 │ 2.64232      3.0  5.57183
  12 │ 2.08189      1.0  8.72759

```
"""
function rand(distribution::DesignDistribution, n::Int = 1)
    RandomDesign(DataFrame(random_design!(values(distribution.factors), n),
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
f1: Distributions.Uniform{Float64}(a=2.0, b=3.0)
f2: Distributions.DiscreteUniform(a=-1, b=5)
f3: Distributions.Uniform{Float64}(a=5.0, b=10.0)

julia> design = rand(design_distribution, 400);

julia> f = @formula 0 ~ f1 + f2 + f3 + f2 ^ 2;

julia> OptimalDesign(design, f, 10)
OptimalDesign
Dimension: (10, 3)
Factors: (f1 = Distributions.Uniform{Float64}(a=2.0, b=3.0), f2 = Distributions.DiscreteUniform(a=-1, b=5), f3 = Distributions.Uniform{Float64}(a=5.0, b=10.0))
Formula: 0 ~ f1 + f2 + f3 + :(f2 ^ 2)
Selected Candidate Rows: [244, 49, 375, 43, 369, 44, 16, 346, 175, 205]
Optimality Criteria: Dict(:D => 2.3940431912232483)
Design Matrix:
10×3 DataFrame
 Row │ f1       f2       f3
     │ Float64  Float64  Float64
─────┼───────────────────────────
   1 │ 2.99329      1.0  5.09246
   2 │ 2.96899      5.0  9.39802
   3 │ 2.96274     -1.0  5.31426
   4 │ 2.00285      2.0  5.40398
   5 │ 2.8491       1.0  9.90621
   6 │ 2.00309     -1.0  9.49394
   7 │ 2.20051      2.0  9.75605
   8 │ 2.06422      5.0  5.1759
   9 │ 2.037       -1.0  9.13114
  10 │ 2.82612      5.0  9.66349

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
