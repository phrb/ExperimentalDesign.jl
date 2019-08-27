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
    design_matrix::DataFrame
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
    matrix::Array{Float64, 2}
    response::Union{Array{Float64, 1}, Missing}
    formula::Union{FormulaTerm, Missing}
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
