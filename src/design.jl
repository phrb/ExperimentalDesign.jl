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
function PlackettBurman(factors::Tuple)
    symbol_factors = Symbol.(factors)
    design = DataFrame(plackettburman(length(factors)))
    dummy_factors = Tuple(Symbol("d" * string(i)) for i = 1:(length(names(design)) -
                                                             length(symbol_factors)))

    design_names = (symbol_factors..., dummy_factors...)
    names!(design, collect(design_names))

    PlackettBurman(design,
                   symbol_factors,
                   dummy_factors,
                   term(:response) ~ sum(term.(design_names)))
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
