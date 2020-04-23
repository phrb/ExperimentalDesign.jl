function Base.show(io::IO, design::PlackettBurman)
    print(io, "$(typeof(design))\n",
          "Dimension: $(size(design.matrix))\n",
          "Factors: $(design.factors)\n",
          "Dummy Factors: $(design.dummy_factors)\n",
          "Formula: $(design.formula)\n",
          "Design Matrix:\n",
          "$(design.matrix)")
end

function Base.show(io::IO, design::FullFactorial)
    if typeof(design.matrix) == Missing
        print(io, "$(typeof(design))\n",
              "Factors: $(design.factors)\n",
              "Formula: $(design.formula)\n",
              "Design Matrix: On demand")
    else
        print(io, "$(typeof(design))\n",
              "Dimension: $(size(design.matrix))\n",
              "Factors: $(design.factors)\n",
              "Formula: $(design.formula)\n",
              "Design Matrix:\n",
              "$(design.matrix)")
    end
end

function Base.show(io::IO, design::RandomDesign)
    print(io, "$(typeof(design))\n",
          "Formula: $(design.formula)\n",
          "Factor Distributions:")

    for factor in keys(design.factors)
        print(io, "\n$(factor): $(design.factors[factor])")
    end
end
