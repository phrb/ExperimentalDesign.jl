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
PlackettBurman(8×7 DataFrames.DataFrame
│ Row │ x1     │ x2     │ x3     │ x4     │ dummy1 │ dummy2 │ dummy3 │
│     │ Symbol │ Symbol │ Symbol │ Symbol │ Symbol │ Symbol │ Symbol │
├─────┼────────┼────────┼────────┼────────┼────────┼────────┼────────┤
│ 1   │ high   │ high   │ high   │ high   │ high   │ high   │ high   │
│ 2   │ low    │ high   │ low    │ high   │ high   │ low    │ low    │
│ 3   │ high   │ low    │ high   │ high   │ low    │ low    │ low    │
│ 4   │ low    │ high   │ high   │ low    │ low    │ low    │ high   │
│ 5   │ high   │ high   │ low    │ low    │ low    │ high   │ low    │
│ 6   │ high   │ low    │ low    │ low    │ high   │ low    │ high   │
│ 7   │ low    │ low    │ low    │ high   │ low    │ high   │ high   │
│ 8   │ low    │ low    │ high   │ low    │ high   │ high   │ low    │, (:x1, :x2, :x3, :x4), (:dummy1, :dummy2, :dummy3), y ~ x1 + x2 + x3 + x4 + dummy1 + dummy2 + dummy3)

```
"""
function PlackettBurman(formula::FormulaTerm)
    symbol_factors = Tuple(r.sym for r in formula.rhs)

    initial_design = plackettburman(length(symbol_factors))
    categorical_design = similar(initial_design, Symbol)

    for i = 1:size(initial_design, 1), j = 1:size(initial_design, 2)
        categorical_design[i, j] = initial_design[i, j] == 1.0 ? :high : :low
    end

    design = DataFrame(categorical_design)

    dummy_factors = Tuple(Symbol("dummy" * string(i)) for i = 1:(length(names(design)) -
                                                                 length(symbol_factors)))

    design_names = (symbol_factors..., dummy_factors...)
    names!(design, collect(design_names))

    PlackettBurman(design,
                   symbol_factors,
                   dummy_factors,
                   term(formula.lhs) ~ sum(term.(design_names)))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> PlackettBurman(4)
PlackettBurman(8×7 DataFrames.DataFrame
│ Row │ factor1 │ factor2 │ factor3 │ factor4 │ dummy1 │ dummy2 │ dummy3 │
│     │ Symbol  │ Symbol  │ Symbol  │ Symbol  │ Symbol │ Symbol │ Symbol │
├─────┼─────────┼─────────┼─────────┼─────────┼────────┼────────┼────────┤
│ 1   │ high    │ high    │ high    │ high    │ high   │ high   │ high   │
│ 2   │ low     │ high    │ low     │ high    │ high   │ low    │ low    │
│ 3   │ high    │ low     │ high    │ high    │ low    │ low    │ low    │
│ 4   │ low     │ high    │ high    │ low     │ low    │ low    │ high   │
│ 5   │ high    │ high    │ low     │ low     │ low    │ high   │ low    │
│ 6   │ high    │ low     │ low     │ low     │ high   │ low    │ high   │
│ 7   │ low     │ low     │ low     │ high    │ low    │ high   │ high   │
│ 8   │ low     │ low     │ high    │ low     │ high   │ high   │ low    │, (:factor1, :factor2, :factor3, :factor4), (:dummy1, :dummy2, :dummy3), response ~ factor1 + factor2 + factor3 + factor4 + dummy1 + dummy2 + dummy3)

```
"""
function PlackettBurman(factors::Int)
    PlackettBurman(term(:response) ~ sum(term.(Symbol("factor" * string(i)) for i = 1:factors)))
end

"""
$(TYPEDEF)

Encapsulates a full  factorial design. Can contain an explicit  design where all
experiments are completely computed and stored, or an iterator of a design, from
which experiments may be extracted as needed.

$(TYPEDFIELDS)
"""
struct FullFactorial <: AbstractFactorialDesign
    matrix::Union{DataFrame, Missing}
    iterator::Base.Iterators.ProductIterator
    factors::NamedTuple
    formula::FormulaTerm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> FullFactorial((A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0]), @formula(y ~ A + B +C), explicit = true)
FullFactorial(12×3 DataFrames.DataFrame
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
│ 12  │ 4   │ b   │ -1.0 │, Base.Iterators.ProductIterator{Tuple{Array{Int64,1},Array{Symbol,1},Array{Float64,1}}}(([1, 2, 4], Symbol[:a, :b], [1.0, -1.0])), (A = [1, 2, 4], B = Symbol[:a, :b], C = [1.0, -1.0]), y ~ A + B + C)

```
"""
function FullFactorial(factors::NamedTuple, formula::FormulaTerm; explicit::Bool = false)
    iterator = fullfactorial(values(factors))
    matrix = missing

    if explicit
        matrix = DataFrame(explicit_fullfactorial(iterator))
        names!(matrix, [r.sym for r in formula.rhs])
    end

    FullFactorial(matrix, iterator, factors, formula)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> FullFactorial((A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0]), explicit = true)
FullFactorial(12×3 DataFrames.DataFrame
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
│ 12  │ 4   │ b   │ -1.0 │, Base.Iterators.ProductIterator{Tuple{Array{Int64,1},Array{Symbol,1},Array{Float64,1}}}(([1, 2, 4], Symbol[:a, :b], [1.0, -1.0])), (A = [1, 2, 4], B = Symbol[:a, :b], C = [1.0, -1.0]), response ~ A + B + C)

```
"""
function FullFactorial(factors::NamedTuple; explicit::Bool = false)
    FullFactorial(factors, term(:response) ~ sum(term.(keys(factors))), explicit = explicit)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> FullFactorial(([1, 2, 4], [:a, :b], [1.0, -1.0]), explicit = true)
FullFactorial(12×3 DataFrames.DataFrame
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
│ 12  │ 4       │ b       │ -1.0    │, Base.Iterators.ProductIterator{Tuple{Array{Int64,1},Array{Symbol,1},Array{Float64,1}}}(([1, 2, 4], Symbol[:a, :b], [1.0, -1.0])), (factor1 = [1, 2, 4], factor2 = Symbol[:a, :b], factor3 = [1.0, -1.0]), response ~ factor1 + factor2 + factor3)

```
"""
function FullFactorial(factors::Tuple; explicit::Bool = false)
    FullFactorial(NamedTuple{Tuple(Symbol("factor" * string(i)) for i = 1:length(factors))}(factors), explicit = explicit)
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

Encapsulates a random design generator.  Receives a `NamedTuple` of factor names
associated  with `Distribution`s  from the  the *Distributions*  package.  After
instantiating  a  `RandomDesign`,  you  must   request  samples  from  it  using
[`rand`](@ref).

"""
struct RandomDesign <: AbstractRandomDesign
    factors::NamedTuple
    formula::FormulaTerm
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> RandomDesign((f1 = Distributions.Uniform(2, 3), f2 = Distributions.DiscreteUniform(-1, 5), f3 = Distributions.Uniform(5, 10)))
RandomDesign((f1 = Distributions.Uniform{Float64}(a=2.0, b=3.0), f2 = Distributions.DiscreteUniform(a=-1, b=5), f3 = Distributions.Uniform{Float64}(a=5.0, b=10.0)), response ~ f1 + f2 + f3)

```
"""
function RandomDesign(factors::NamedTuple)
    RandomDesign(factors, term(:response) ~ sum(term.(keys(factors))))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> RandomDesign((Distributions.Uniform(2, 3), Distributions.DiscreteUniform(-1, 5), Distributions.Uniform(5, 10)))
RandomDesign((factor1 = Distributions.Uniform{Float64}(a=2.0, b=3.0), factor2 = Distributions.DiscreteUniform(a=-1, b=5), factor3 = Distributions.Uniform{Float64}(a=5.0, b=10.0)), response ~ factor1 + factor2 + factor3)

```
"""
function RandomDesign(factors::Tuple)
    RandomDesign(NamedTuple{Tuple(Symbol("factor" * string(i)) for i = 1:length(factors))}(factors))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> rand(RandomDesign((f1 = Distributions.Uniform(2, 3), f2 = Distributions.DiscreteUniform(-1, 5), f3 = Distributions.Uniform(5, 10))), 12)
12×3 DataFrames.DataFrame
│ Row │ f1      │ f2   │ f3      │
│     │ Real    │ Real │ Real    │
├─────┼─────────┼──────┼─────────┤
│ 1   │ 2.04922 │ 5    │ 8.85741 │
│ 2   │ 2.25659 │ -1   │ 6.57617 │
│ 3   │ 2.86526 │ 4    │ 5.41422 │
│ 4   │ 2.81201 │ -1   │ 5.40944 │
│ 5   │ 2.92412 │ 1    │ 6.22902 │
│ 6   │ 2.20051 │ 5    │ 8.80749 │
│ 7   │ 2.69459 │ 3    │ 6.34626 │
│ 8   │ 2.28902 │ 2    │ 8.72759 │
│ 9   │ 2.95173 │ 0    │ 5.42517 │
│ 10  │ 2.35163 │ 4    │ 5.89413 │
│ 11  │ 2.97614 │ 5    │ 9.99911 │
│ 12  │ 2.4577  │ 0    │ 8.303   │

```
"""
function rand(design::RandomDesign, n::Int = 1)
    DataFrame(random_design(values(design.factors), n), collect(keys(design.factors)))
end


"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct OptimalDesign <: AbstractOptimalDesign
    matrix::DataFrame
    candidates::AbstractDesign
    criteria::Dict{Symbol, Float64}
    formula::FormulaTerm
end
