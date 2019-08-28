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
"""
function PlackettBurman(factors::Int)
    PlackettBurman(term(:response) ~ sum(term.(Symbol("factor" * string(i)) for i = 1:factors)))
end

"""
$(TYPEDEF)

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
    matrix::Array{Float64, 2}
    response::Union{Array{Float64, 1}, Missing}
    formula::Union{FormulaTerm, Missing}
    aliases::Union{FormulaTerm, Missing}
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct Optimal <: AbstractOptimalDesign
    matrix::Array{Float64, 2}
    response::Union{Array{Float64, 1}, Missing}
    formula::Union{FormulaTerm, Missing}
    optimality::Union{NamedTuple, Missing}
end
