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
abstract type AbstractResponseSurfaceDesign <: AbstractDesign end

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

Box-Behnken design for response surface methodology.

$(TYPEDFIELDS)
"""
struct BoxBehnken <: AbstractResponseSurfaceDesign
    matrix::DataFrame
    factors::Tuple
    formula::FormulaTerm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> BoxBehnken(@formula( y ~ -1 + factor1 & factor1 + factor2 & factor2 + factor3 & factor3 + factor1 & factor2 + factor1 & factor3 + factor2 & factor3 + factor1 + factor2 + factor3), center=3)
BoxBehnken
Dimension: (15, 3)
Factors: (:factor1, :factor2, :factor3)
Formula: y ~ -1 + factor1 + factor2 + factor3 + factor1 & factor1 + factor2 & factor2 + factor3 & factor3 + factor1 & factor2 + factor1 & factor3 + factor2 & factor3
Design Matrix:
15×3 DataFrame
 Row │ factor1  factor2  factor3
     │ Float64  Float64  Float64
─────┼───────────────────────────
   1 │    -1.0     -1.0      0.0
   2 │     1.0     -1.0      0.0
   3 │    -1.0      1.0      0.0
   4 │     1.0      1.0      0.0
   5 │    -1.0      0.0     -1.0
   6 │     1.0      0.0     -1.0
   7 │    -1.0      0.0      1.0
   8 │     1.0      0.0      1.0
   9 │     0.0     -1.0     -1.0
  10 │     0.0      1.0     -1.0
  11 │     0.0     -1.0      1.0
  12 │     0.0      1.0      1.0
  13 │     0.0      0.0      0.0
  14 │     0.0      0.0      0.0
  15 │     0.0      0.0      0.0

```
"""
function BoxBehnken(formula::FormulaTerm; center::Int=1)
    symbol_factors_array = []
    # TODO can only handle interaction terms not functional terms
    for r in formula.rhs
        if r isa InteractionTerm
            for s in r.terms
                append!(symbol_factors_array, [s.sym])
            end
        elseif r isa Term
            append!(symbol_factors_array, [r.sym])
        end
    end
    symbol_factors = Tuple(unique!(symbol_factors_array))

    initial_design = boxbehnken(length(symbol_factors_array), center)
    design = DataFrame(initial_design, :auto)

    rename!(design, collect(symbol_factors))

    BoxBehnken(design, symbol_factors, formula)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> BoxBehnken(3, center=3)
BoxBehnken
Dimension: (15, 3)
Factors: (:factor1, :factor2, :factor3)
Formula: y ~ -1 + factor1 & factor1 + factor2 & factor2 + factor3 & factor3 + factor1 & factor2 + factor1 & factor3 + factor2 & factor3 + factor1 + factor2 + factor3
Design Matrix:
15×3 DataFrame
 Row │ factor1  factor2  factor3
     │ Float64  Float64  Float64
─────┼───────────────────────────
   1 │    -1.0     -1.0      0.0
   2 │     1.0     -1.0      0.0
   3 │    -1.0      1.0      0.0
   4 │     1.0      1.0      0.0
   5 │    -1.0      0.0     -1.0
   6 │     1.0      0.0     -1.0
   7 │    -1.0      0.0      1.0
   8 │     1.0      0.0      1.0
   9 │     0.0     -1.0     -1.0
  10 │     0.0      1.0     -1.0
  11 │     0.0     -1.0      1.0
  12 │     0.0      1.0      1.0
  13 │     0.0      0.0      0.0
  14 │     0.0      0.0      0.0
  15 │     0.0      0.0      0.0

```
"""
function BoxBehnken(factors::Int; center::Int=1)
    factors_tuple = Tuple(Symbol("factor" * string(i)) for i = 1:factors)
    q_terms = sum(InteractionTerm(term.((only(x),only(x)))) for x ∈ combinations(factors_tuple, 1))
    i_terms = sum((&)(term.(x)...) for x ∈ combinations(factors_tuple, 2))
    l_terms = sum(term(only(x)) for x ∈ combinations(factors_tuple, 1))
    formula = term(:y) ~ term(-1) + q_terms + i_terms + l_terms
    BoxBehnken(formula, center=center)
end

"""
$(TYPEDEF)

Encapsulates a central composite design.

$(TYPEDFIELDS)
"""
struct CentralComposite <: AbstractResponseSurfaceDesign
    matrix::DataFrame
    factors::Tuple
    formula::FormulaTerm
    alpha::Symbol
    face::Symbol
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> CentralComposite(@formula(y ~ -1 + factor1 & factor1 + factor2 & factor2 + factor3 & factor3 + factor1 & factor2 + factor1 & factor3 + factor2 & factor3 + factor1 + factor2 + factor3))
CentralComposite
Dimension: (22, 3)
Factors: (:factor1, :factor2, :factor3)
Formula: y ~ -1 + factor1 + factor2 + factor3 + factor1 & factor1 + factor2 & factor2 + factor3 & factor3 + factor1 & factor2 + factor1 & factor3 + factor2 & factor3
Alpha: orthogonal
Face: circumscribed
Design Matrix:
22×3 DataFrame
 Row │ factor1   factor2   factor3
     │ Float64   Float64   Float64
─────┼──────────────────────────────
   1 │ -1.0      -1.0      -1.0
   2 │  1.0      -1.0      -1.0
   3 │ -1.0       1.0      -1.0
   4 │  1.0       1.0      -1.0
   5 │ -1.0      -1.0       1.0
   6 │  1.0      -1.0       1.0
   7 │ -1.0       1.0       1.0
   8 │  1.0       1.0       1.0
  ⋮  │    ⋮         ⋮         ⋮
  16 │  0.0       1.82574   0.0
  17 │  0.0       0.0      -1.82574
  18 │  0.0       0.0       1.82574
  19 │  0.0       0.0       0.0
  20 │  0.0       0.0       0.0
  21 │  0.0       0.0       0.0
  22 │  0.0       0.0       0.0
                      7 rows omitted

```
"""
function CentralComposite(formula::FormulaTerm; center::Array{Int}=[4, 4], alpha::Symbol=:orthogonal, face::Symbol=:circumscribed)
    symbol_factors_array = []
    # TODO can only handle interaction terms not functional terms
    for r in formula.rhs
        if r isa InteractionTerm
            for s in r.terms
                append!(symbol_factors_array, [s.sym])
            end
        elseif r isa Term
            append!(symbol_factors_array, [r.sym])
        end
    end
    symbol_factors = Tuple(unique!(symbol_factors_array))

    initial_design = ccdesign(length(symbol_factors_array), center, alpha, face)
    design = DataFrame(initial_design, :auto)
    rename!(design, collect(symbol_factors))

    CentralComposite(design, symbol_factors, formula, alpha, face)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> CentralComposite(3)
CentralComposite
Dimension: (22, 3)
Factors: (:factor1, :factor2, :factor3)
Formula: y ~ -1 + factor1 & factor1 + factor2 & factor2 + factor3 & factor3 + factor1 & factor2 + factor1 & factor3 + factor2 & factor3 + factor1 + factor2 + factor3
Alpha: orthogonal
Face: circumscribed
Design Matrix:
22×3 DataFrame
 Row │ factor1   factor2   factor3
     │ Float64   Float64   Float64
─────┼──────────────────────────────
   1 │ -1.0      -1.0      -1.0
   2 │  1.0      -1.0      -1.0
   3 │ -1.0       1.0      -1.0
   4 │  1.0       1.0      -1.0
   5 │ -1.0      -1.0       1.0
   6 │  1.0      -1.0       1.0
   7 │ -1.0       1.0       1.0
   8 │  1.0       1.0       1.0
  ⋮  │    ⋮         ⋮         ⋮
  16 │  0.0       1.82574   0.0
  17 │  0.0       0.0      -1.82574
  18 │  0.0       0.0       1.82574
  19 │  0.0       0.0       0.0
  20 │  0.0       0.0       0.0
  21 │  0.0       0.0       0.0
  22 │  0.0       0.0       0.0
                      7 rows omitted

```
"""
function CentralComposite(factors::Int; center::Array{Int}=[4, 4], alpha::Symbol=:orthogonal, face::Symbol=:circumscribed)
    factors_tuple = Tuple(Symbol("factor" * string(i)) for i = 1:factors)
    q_terms = sum(InteractionTerm(term.((only(x),only(x)))) for x ∈ combinations(factors_tuple, 1))
    i_terms = sum((&)(term.(x)...) for x ∈ combinations(factors_tuple, 2))
    l_terms = sum(term(only(x)) for x ∈ combinations(factors_tuple, 1))
    formula = term(:y) ~ term(-1) + q_terms + i_terms + l_terms
    CentralComposite(formula, center=center,  alpha=alpha, face=face)
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
struct FractionalFactorial2Level <: AbstractFactorialDesign
    matrix::DataFrame
    factors::NamedTuple
    formula::FormulaTerm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> FractionalFactorial2Level(@formula(y ~ a + b + a&b + c + a&c+ b&c + a&b&c))
FractionalFactorial2Level
Dimension: (8, 7)
Factors: (a = [-1, 1], b = [-1, 1], c = [-1, 1])
Formula: y ~ a + b + c + a & b + a & c + b & c + a & b & c
Design Matrix:
8×7 DataFrame
 Row │ a      b      c      a_b    a_c    b_c    a_b_c
     │ Int64  Int64  Int64  Int64  Int64  Int64  Int64
─────┼─────────────────────────────────────────────────
   1 │    -1     -1     -1      1      1      1     -1
   2 │     1     -1     -1     -1     -1      1      1
   3 │    -1      1     -1     -1      1     -1      1
   4 │     1      1     -1      1     -1     -1     -1
   5 │    -1     -1      1      1     -1     -1      1
   6 │     1     -1      1     -1      1     -1     -1
   7 │    -1      1      1     -1     -1      1     -1
   8 │     1      1      1      1      1      1      1

```
"""
function FractionalFactorial2Level(formula::FormulaTerm)
    symbol_main_factors_array = []
    symbol_interaction_factors_array = []
    # TODO can only handle simple interaction terms but not functional terms
    for r in formula.rhs
        if r isa InteractionTerm
            append!(symbol_interaction_factors_array, [r])
            for s in r.terms
                append!(symbol_main_factors_array, [s.sym])
            end
        elseif r isa Term
            append!(symbol_main_factors_array, [r.sym])
        end
    end
    symbol_factors = Tuple(unique!(symbol_main_factors_array))
    levels = [[-1,1] for i in symbol_factors]
    A_ff = FullFactorial(NamedTuple{symbol_factors}(levels))
    for (index, interactionterm) in enumerate(symbol_interaction_factors_array)
        interaction_factors = [term.sym for term in interactionterm.terms]
        switch = false
        fac = A_ff.matrix[:, interaction_factors[1]]
        for factor in interaction_factors
            if switch == false
                switch = true
            else
                fac =  fac .* A_ff.matrix[!, factor]
            end
        end
        A_ff.matrix[!, Symbol(join(interaction_factors, "_"))] = fac
    end
    FractionalFactorial2Level(A_ff.matrix, NamedTuple{symbol_factors}(levels), formula)
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
factor1: DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor2: DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor3: DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor4: DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor5: DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])
factor6: DiscreteNonParametric{Int64, Float64, Vector{Int64}, Vector{Float64}}(support=[-1, 1], p=[0.5, 0.5])

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
"""
struct RandomLHCDesign <: AbstractRandomDesign
    matrix::DataFrame
    factors::Tuple
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> RandomLHCDesign(4,3)
RandomLHCDesign
Dimension: (4, 3)
Factors: (:factor1, :factor2, :factor3)
Design Matrix:
4×3 DataFrame
 Row │ factor1  factor2  factor3
     │ Int64    Int64    Int64
─────┼───────────────────────────
   1 │       1        4        4
   2 │       2        1        2
   3 │       3        3        3
   4 │       4        2        1

```
"""
function RandomLHCDesign(n::Int, d::Int)
    factors = Tuple(Symbol("factor" * string(i)) for i = 1:d)
    design = DataFrame(randomLHC(n,d), :auto)
    rename!(design, collect(factors))
    RandomLHCDesign(design, factors)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> RandomLHCDesign(4,Tuple((:A,:B,:C)))
RandomLHCDesign
Dimension: (4, 3)
Factors: (:A, :B, :C)
Design Matrix:
4×3 DataFrame
 Row │ A      B      C
     │ Int64  Int64  Int64
─────┼─────────────────────
   1 │     1      4      4
   2 │     2      1      2
   3 │     3      3      3
   4 │     4      2      1
```
"""
function RandomLHCDesign(n::Int, factors::Tuple)
    design = DataFrame(randomLHC(n,length(factors)), :auto)
    rename!(design, collect(factors))
    RandomLHCDesign(design, factors)
end


"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct OptimLHCDesign <: AbstractRandomDesign
    matrix::DataFrame
    factors::Tuple
    fitness::Vector{Float64}
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> OptimLHCDesign(5,3,4)
OptimLHCDesign
Dimension: (5, 3)
Factors: (:factor1, :factor2, :factor3)
Fitness: [1.3760238272524201, 1.3760238272524201, 1.3760238272524201, 1.3760238272524201, 1.3760238272524201]
Design Matrix:
5×3 DataFrame
 Row │ factor1  factor2  factor3
     │ Int64    Int64    Int64
─────┼───────────────────────────
   1 │       3        5        2
   2 │       5        1        3
   3 │       4        4        5
   4 │       2        2        1
   5 │       1        3        4
```
"""
function OptimLHCDesign(n::Int, d::Int, gens)
    factors = Tuple(Symbol("factor" * string(i)) for i = 1:d)
    design_matrix, fitness = LHCoptim(n, d, gens)
    design = DataFrame(design_matrix, :auto)
    rename!(design, collect(factors))
    OptimLHCDesign(design, factors, fitness)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> OptimLHCDesign(5,Tuple((:A,:B,:C)),4)
OptimLHCDesign
Dimension: (5, 3)
Factors: (:A, :B, :C)
Fitness: [1.3760238272524201, 1.3760238272524201, 1.3760238272524201, 1.3760238272524201, 1.3760238272524201]
Design Matrix:
5×3 DataFrame
 Row │ A      B      C
     │ Int64  Int64  Int64
─────┼─────────────────────
   1 │     3      5      2
   2 │     5      1      3
   3 │     4      4      5
   4 │     2      2      1
   5 │     1      3      4

```
"""
function OptimLHCDesign(n::Int, factors::Tuple, gens)
    design_matrix, fitness = LHCoptim(n, length(factors), gens)
    design = DataFrame(design_matrix, :auto)
    rename!(design, collect(factors))
    OptimLHCDesign(design, factors, fitness)
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
