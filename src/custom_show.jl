function Base.show(io::IO, design::PlackettBurman)
    print(io, "$(typeof(design))\n",
          "Dimension: $(size(design.matrix))\n",
          "Factors: $(design.factors)\n",
          "Dummy Factors: $(design.dummy_factors)\n",
          "Formula: $(design.formula)\n",
          "Design Matrix:\n")
    show(io, design.matrix)
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
              "Design Matrix:\n")
        show(io, design.matrix)
    end
end

function Base.show(io::IO, design::FractionalFactorial2Level)
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
              "Design Matrix:\n")
        show(io, design.matrix)
    end
end

function Base.show(io::IO, design::DesignDistribution)
    print(io, "$(typeof(design))\n",
          "Formula: $(design.formula)\n",
          "Factor Distributions:")

    for factor in keys(design.factors)
        print(io, "\n$(factor): $(design.factors[factor])")
    end
end

function Base.show(io::IO, design::RandomDesign)
    print(io, "$(typeof(design))\n",
          "Dimension: $(size(design.matrix))\n",
          "Factors: $(design.factors)\n",
          "Formula: $(design.formula)\n",
          "Design Matrix:\n")
    show(io, design.matrix)
end

function Base.show(io::IO, design::RandomLHCDesign)
    print(io, "$(typeof(design))\n",
          "Dimension: $(size(design.matrix))\n",
          "Factors: $(design.factors)\n",
          "Design Matrix:\n")
    show(io, design.matrix)
end

function Base.show(io::IO, design::OptimLHCDesign)
    print(io, "$(typeof(design))\n",
          "Dimension: $(size(design.matrix))\n",
          "Factors: $(design.factors)\n",
          "Fitness: $(design.fitness)\n",
          "Design Matrix:\n")
    show(io, design.matrix)
end

function Base.show(io::IO, design::OptimalDesign)
    print(io, "$(typeof(design))\n",
          "Dimension: $(size(design.matrix))\n",
          "Factors: $(design.factors)\n",
          "Formula: $(design.formula)\n",
          "Selected Candidate Rows: $(design.selected_experiments)\n",
          "Optimality Criteria: $(design.criteria)\n",
          "Design Matrix:\n")
    show(io, design.matrix)
end
