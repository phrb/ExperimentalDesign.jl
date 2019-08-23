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
struct PlackettBurmanDesign <: AbstractScreeningDesign
    matrix::Array{Float64, 2}
    formula::Union{FormulaTerm, Missing}
    aliases::Union{FormulaTerm, Missing}
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct FactorialDesign <: AbstractFactorialDesign
    matrix::Array{Float64, 2}
    formula::Union{FormulaTerm, Missing}
    resolution::Union{Int, Missing}
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct OptimalDesign <: AbstractOptimalDesign
    matrix::Array{Float64, 2}
    formula::Union{FormulaTerm, Missing}
    optimality::Union{NamedTuple, Missing}
end
