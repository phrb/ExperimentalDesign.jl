var documenterSearchIndex = {"docs":
[{"location":"lib/internals/#API-Documentation-1","page":"Internals","title":"API Documentation","text":"","category":"section"},{"location":"lib/internals/#","page":"Internals","title":"Internals","text":"Documentation for ExperimentalDesign.jl's API.","category":"page"},{"location":"lib/internals/#Contents-1","page":"Internals","title":"Contents","text":"","category":"section"},{"location":"lib/internals/#","page":"Internals","title":"Internals","text":"Pages = [\"internals.md\"]","category":"page"},{"location":"lib/internals/#Index-1","page":"Internals","title":"Index","text":"","category":"section"},{"location":"lib/internals/#","page":"Internals","title":"Internals","text":"Pages = [\"internals.md\"]","category":"page"},{"location":"lib/internals/#API-1","page":"Internals","title":"API","text":"","category":"section"},{"location":"lib/internals/#","page":"Internals","title":"Internals","text":"Modules = [ExperimentalDesign]\nPublic = false","category":"page"},{"location":"lib/internals/#ExperimentalDesign.AbstractRandomDesign","page":"Internals","title":"ExperimentalDesign.AbstractRandomDesign","text":"abstract type AbstractRandomDesign <: AbstractDesign\n\n\n\n\n\n","category":"type"},{"location":"lib/internals/#ExperimentalDesign.IterableFullFactorial","page":"Internals","title":"ExperimentalDesign.IterableFullFactorial","text":"struct IterableFullFactorial <: AbstractFactorialDesign\n\nmatrix::DataFrames.DataFrame\nfactors::NamedTuple\nformula::StatsModels.FormulaTerm\n\n\n\n\n\n","category":"type"},{"location":"lib/internals/#ExperimentalDesign.RandomDesign","page":"Internals","title":"ExperimentalDesign.RandomDesign","text":"struct RandomDesign <: ExperimentalDesign.AbstractRandomDesign\n\nmatrix::DataFrames.DataFrame\nfactors::NamedTuple\nformula::StatsModels.FormulaTerm\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#API-Documentation-1","page":"Public","title":"API Documentation","text":"","category":"section"},{"location":"lib/public/#","page":"Public","title":"Public","text":"Documentation for ExperimentalDesign.jl's API.","category":"page"},{"location":"lib/public/#Contents-1","page":"Public","title":"Contents","text":"","category":"section"},{"location":"lib/public/#","page":"Public","title":"Public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"lib/public/#Index-1","page":"Public","title":"Index","text":"","category":"section"},{"location":"lib/public/#","page":"Public","title":"Public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"lib/public/#API-1","page":"Public","title":"API","text":"","category":"section"},{"location":"lib/public/#","page":"Public","title":"Public","text":"Modules = [ExperimentalDesign]\nPrivate = false","category":"page"},{"location":"lib/public/#ExperimentalDesign.AbstractDesign","page":"Public","title":"ExperimentalDesign.AbstractDesign","text":"abstract type AbstractDesign\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.AbstractFactorialDesign","page":"Public","title":"ExperimentalDesign.AbstractFactorialDesign","text":"abstract type AbstractFactorialDesign <: AbstractDesign\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.AbstractOptimalDesign","page":"Public","title":"ExperimentalDesign.AbstractOptimalDesign","text":"abstract type AbstractOptimalDesign <: AbstractDesign\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.AbstractScreeningDesign","page":"Public","title":"ExperimentalDesign.AbstractScreeningDesign","text":"abstract type AbstractScreeningDesign <: AbstractDesign\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.CategoricalFactor","page":"Public","title":"ExperimentalDesign.CategoricalFactor","text":"struct CategoricalFactor <: Distributions.Distribution{Distributions.Univariate,Distributions.Discrete}\n\nvalues::Array{Any,1}\ndistribution::Distributions.DiscreteUniform\n\nA simple wrapper for a DiscreteUniform distribution over non-numerical arrays.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.CategoricalFactor-Tuple{Array{Any,1}}","page":"Public","title":"ExperimentalDesign.CategoricalFactor","text":"CategoricalFactor(values::Array{Any,1}) -> CategoricalFactor\n\n\njulia> a = CategoricalFactor([:a, :b, 2, 1.0])\nCategoricalFactor(\nvalues: Any[:a, :b, 2, 1.0]\ndistribution: Distributions.DiscreteUniform(a=1, b=4)\n)\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.DesignDistribution","page":"Public","title":"ExperimentalDesign.DesignDistribution","text":"struct DesignDistribution\n\nfactors::NamedTuple\nformula::StatsModels.FormulaTerm\n\nEncapsulates one or more Distributions, generating random designs.  Receives a NamedTuple  of  factor names  associated  with  Distributions from  the  the Distributions package.   After instantiating a DesignDistribution,  you must request   samples   from   it    using   rand,   which   generates   a RandomDesign\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.DesignDistribution-Tuple{Array}","page":"Public","title":"ExperimentalDesign.DesignDistribution","text":"DesignDistribution(factors::Array) -> DesignDistribution\n\n\njulia> DesignDistribution([Uniform(2, 3), DiscreteUniform(-1, 5), Uniform(5, 10)])\nDesignDistribution\nFormula: 0 ~ factor1 + factor2 + factor3\nFactor Distributions:\nfactor1: Distributions.Uniform{Float64}(a=2.0, b=3.0)\nfactor2: Distributions.DiscreteUniform(a=-1, b=5)\nfactor3: Distributions.Uniform{Float64}(a=5.0, b=10.0)\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.DesignDistribution-Tuple{NamedTuple}","page":"Public","title":"ExperimentalDesign.DesignDistribution","text":"DesignDistribution(factors::NamedTuple) -> DesignDistribution\n\n\njulia> DesignDistribution((f1 = Uniform(2, 3), f2 = DiscreteUniform(-1, 5), f3 = Uniform(5, 10)))\nDesignDistribution\nFormula: 0 ~ f1 + f2 + f3\nFactor Distributions:\nf1: Distributions.Uniform{Float64}(a=2.0, b=3.0)\nf2: Distributions.DiscreteUniform(a=-1, b=5)\nf3: Distributions.Uniform{Float64}(a=5.0, b=10.0)\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.DesignDistribution-Tuple{Tuple}","page":"Public","title":"ExperimentalDesign.DesignDistribution","text":"DesignDistribution(factors::Tuple) -> DesignDistribution\n\n\njulia> DesignDistribution((Uniform(2, 3), DiscreteUniform(-1, 5), Uniform(5, 10)))\nDesignDistribution\nFormula: 0 ~ factor1 + factor2 + factor3\nFactor Distributions:\nfactor1: Distributions.Uniform{Float64}(a=2.0, b=3.0)\nfactor2: Distributions.DiscreteUniform(a=-1, b=5)\nfactor3: Distributions.Uniform{Float64}(a=5.0, b=10.0)\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.DesignDistribution-Union{Tuple{D}, Tuple{D,Int64}} where D<:Distributions.Distribution","page":"Public","title":"ExperimentalDesign.DesignDistribution","text":"DesignDistribution(distribution::D<:Distributions.Distribution, n::Int64) -> DesignDistribution\n\n\njulia> DesignDistribution(DiscreteNonParametric([-1, 1], [0.5, 0.5]), 6)\nDesignDistribution\nFormula: 0 ~ factor1 + factor2 + factor3 + factor4 + factor5 + factor6\nFactor Distributions:\nfactor1: Distributions.DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])\nfactor2: Distributions.DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])\nfactor3: Distributions.DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])\nfactor4: Distributions.DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])\nfactor5: Distributions.DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])\nfactor6: Distributions.DiscreteNonParametric{Int64,Float64,Array{Int64,1},Array{Float64,1}}(support=[-1, 1], p=[0.5, 0.5])\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.FractionalFactorial","page":"Public","title":"ExperimentalDesign.FractionalFactorial","text":"struct FractionalFactorial <: AbstractFactorialDesign\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.FullFactorial","page":"Public","title":"ExperimentalDesign.FullFactorial","text":"struct FullFactorial <: AbstractFactorialDesign\n\nEncapsulates  an explicit  full factorial  design, where  all experiments  are completely computed  and stored.  Use IterableFullFactorial  to obtain an implicit factorial  design with arbitrary access to the  elements of a full factorial design.\n\nmatrix::DataFrames.DataFrame\nfactors::NamedTuple\nformula::StatsModels.FormulaTerm\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.FullFactorial-Tuple{Array}","page":"Public","title":"ExperimentalDesign.FullFactorial","text":"FullFactorial(factors::Array) -> FullFactorial\n\n\njulia> FullFactorial(fill([-1, 1], 3))\nFullFactorial\nDimension: (8, 3)\nFactors: (factor1 = [-1, 1], factor2 = [-1, 1], factor3 = [-1, 1])\nFormula: 0 ~ factor1 + factor2 + factor3\nDesign Matrix:\n8×3 DataFrame\n│ Row │ factor1 │ factor2 │ factor3 │\n│     │ Int64   │ Int64   │ Int64   │\n├─────┼─────────┼─────────┼─────────┤\n│ 1   │ -1      │ -1      │ -1      │\n│ 2   │ 1       │ -1      │ -1      │\n│ 3   │ -1      │ 1       │ -1      │\n│ 4   │ 1       │ 1       │ -1      │\n│ 5   │ -1      │ -1      │ 1       │\n│ 6   │ 1       │ -1      │ 1       │\n│ 7   │ -1      │ 1       │ 1       │\n│ 8   │ 1       │ 1       │ 1       │\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.FullFactorial-Tuple{NamedTuple,StatsModels.FormulaTerm}","page":"Public","title":"ExperimentalDesign.FullFactorial","text":"FullFactorial(factors::NamedTuple, formula::StatsModels.FormulaTerm) -> FullFactorial\n\n\njulia> FullFactorial((A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0]), @formula(y ~ A + B +C))\nFullFactorial\nDimension: (12, 3)\nFactors: (A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0])\nFormula: y ~ A + B + C\nDesign Matrix:\n12×3 DataFrame\n│ Row │ A   │ B   │ C    │\n│     │ Any │ Any │ Any  │\n├─────┼─────┼─────┼──────┤\n│ 1   │ 1   │ a   │ 1.0  │\n│ 2   │ 2   │ a   │ 1.0  │\n│ 3   │ 4   │ a   │ 1.0  │\n│ 4   │ 1   │ b   │ 1.0  │\n│ 5   │ 2   │ b   │ 1.0  │\n│ 6   │ 4   │ b   │ 1.0  │\n│ 7   │ 1   │ a   │ -1.0 │\n│ 8   │ 2   │ a   │ -1.0 │\n│ 9   │ 4   │ a   │ -1.0 │\n│ 10  │ 1   │ b   │ -1.0 │\n│ 11  │ 2   │ b   │ -1.0 │\n│ 12  │ 4   │ b   │ -1.0 │\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.FullFactorial-Tuple{NamedTuple}","page":"Public","title":"ExperimentalDesign.FullFactorial","text":"FullFactorial(factors::NamedTuple) -> FullFactorial\n\n\njulia> FullFactorial((A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0]))\nFullFactorial\nDimension: (12, 3)\nFactors: (A = [1, 2, 4], B = [:a, :b], C = [1.0, -1.0])\nFormula: 0 ~ A + B + C\nDesign Matrix:\n12×3 DataFrame\n│ Row │ A   │ B   │ C    │\n│     │ Any │ Any │ Any  │\n├─────┼─────┼─────┼──────┤\n│ 1   │ 1   │ a   │ 1.0  │\n│ 2   │ 2   │ a   │ 1.0  │\n│ 3   │ 4   │ a   │ 1.0  │\n│ 4   │ 1   │ b   │ 1.0  │\n│ 5   │ 2   │ b   │ 1.0  │\n│ 6   │ 4   │ b   │ 1.0  │\n│ 7   │ 1   │ a   │ -1.0 │\n│ 8   │ 2   │ a   │ -1.0 │\n│ 9   │ 4   │ a   │ -1.0 │\n│ 10  │ 1   │ b   │ -1.0 │\n│ 11  │ 2   │ b   │ -1.0 │\n│ 12  │ 4   │ b   │ -1.0 │\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.FullFactorial-Tuple{Tuple}","page":"Public","title":"ExperimentalDesign.FullFactorial","text":"FullFactorial(factors::Tuple) -> FullFactorial\n\n\njulia> FullFactorial(([1, 2, 4], [:a, :b], [1.0, -1.0]))\nFullFactorial\nDimension: (12, 3)\nFactors: (factor1 = [1, 2, 4], factor2 = [:a, :b], factor3 = [1.0, -1.0])\nFormula: 0 ~ factor1 + factor2 + factor3\nDesign Matrix:\n12×3 DataFrame\n│ Row │ factor1 │ factor2 │ factor3 │\n│     │ Any     │ Any     │ Any     │\n├─────┼─────────┼─────────┼─────────┤\n│ 1   │ 1       │ a       │ 1.0     │\n│ 2   │ 2       │ a       │ 1.0     │\n│ 3   │ 4       │ a       │ 1.0     │\n│ 4   │ 1       │ b       │ 1.0     │\n│ 5   │ 2       │ b       │ 1.0     │\n│ 6   │ 4       │ b       │ 1.0     │\n│ 7   │ 1       │ a       │ -1.0    │\n│ 8   │ 2       │ a       │ -1.0    │\n│ 9   │ 4       │ a       │ -1.0    │\n│ 10  │ 1       │ b       │ -1.0    │\n│ 11  │ 2       │ b       │ -1.0    │\n│ 12  │ 4       │ b       │ -1.0    │\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.OptimalDesign","page":"Public","title":"ExperimentalDesign.OptimalDesign","text":"struct OptimalDesign <: AbstractOptimalDesign\n\nmatrix::DataFrames.DataFrame\nfactors::NamedTuple\nformula::StatsModels.FormulaTerm\nselected_experiments::Array{Int64,N} where N\ncriteria::Dict{Symbol,Float64}\n\nContains a  set of candidate experiments,  a target optimality criterion,  and a model prior  formula.  A set  of experiments  that maximizes a  given optimality criterion can selected from a candidate set using the kl_exchange algorithm.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.OptimalDesign-Tuple{AbstractDesign,StatsModels.FormulaTerm,Int64}","page":"Public","title":"ExperimentalDesign.OptimalDesign","text":"OptimalDesign(candidates::AbstractDesign, formula::StatsModels.FormulaTerm, experiments::Int64; tolerance, seed_design_size, max_iterations, design_k, candidates_l) -> OptimalDesign\n\n\njulia> design_distribution = DesignDistribution((f1 = Uniform(2, 3), f2 = DiscreteUniform(-1, 5), f3 = Uniform(5, 10)))\nDesignDistribution\nFormula: 0 ~ f1 + f2 + f3\nFactor Distributions:\nf1: Distributions.Uniform{Float64}(a=2.0, b=3.0)\nf2: Distributions.DiscreteUniform(a=-1, b=5)\nf3: Distributions.Uniform{Float64}(a=5.0, b=10.0)\n\njulia> design = rand(design_distribution, 400);\n\njulia> f = @formula 0 ~ f1 + f2 + f3 + f2 ^ 2;\n\njulia> OptimalDesign(design, f, 10)\nOptimalDesign\nDimension: (10, 3)\nFactors: (f1 = Distributions.Uniform{Float64}(a=2.0, b=3.0), f2 = Distributions.DiscreteUniform(a=-1, b=5), f3 = Distributions.Uniform{Float64}(a=5.0, b=10.0))\nFormula: 0 ~ f1 + f2 + f3 + :(f2 ^ 2)\nSelected Candidate Rows: [222, 103, 384, 140, 139, 156, 63, 184, 169, 54]\nOptimality Criteria: Dict(:D => 2.418045723405036)\nDesign Matrix:\n10×3 DataFrame\n│ Row │ f1      │ f2      │ f3      │\n│     │ Float64 │ Float64 │ Float64 │\n├─────┼─────────┼─────────┼─────────┤\n│ 1   │ 2.02407 │ -1.0    │ 5.07959 │\n│ 2   │ 2.0203  │ 5.0     │ 9.50193 │\n│ 3   │ 2.80676 │ -1.0    │ 9.87726 │\n│ 4   │ 2.9996  │ 2.0     │ 5.03174 │\n│ 5   │ 2.14846 │ 2.0     │ 9.72596 │\n│ 6   │ 2.84392 │ 5.0     │ 5.01213 │\n│ 7   │ 2.04046 │ 2.0     │ 5.00331 │\n│ 8   │ 2.93432 │ 1.0     │ 9.71405 │\n│ 9   │ 2.79858 │ -1.0    │ 5.35685 │\n│ 10  │ 2.97635 │ 5.0     │ 8.10695 │\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.PlackettBurman","page":"Public","title":"ExperimentalDesign.PlackettBurman","text":"struct PlackettBurman <: AbstractScreeningDesign\n\nEncapsulates  a  Plackett-Burman  screening  design  constructed  using  Paley's method.  Factor  levels are  encoded as  :high and  :low symbols,  and extra dummy variables, possibly  aliasing interactions between factors,  will be added to pad a design.\n\nmatrix::DataFrames.DataFrame\nfactors::Tuple\ndummy_factors::Tuple\nformula::StatsModels.FormulaTerm\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#ExperimentalDesign.PlackettBurman-Tuple{Int64}","page":"Public","title":"ExperimentalDesign.PlackettBurman","text":"PlackettBurman(factors::Int64) -> PlackettBurman\n\n\njulia> PlackettBurman(4)\nPlackettBurman\nDimension: (8, 7)\nFactors: (:factor1, :factor2, :factor3, :factor4)\nDummy Factors: (:dummy1, :dummy2, :dummy3)\nFormula: 0 ~ -1 + factor1 + factor2 + factor3 + factor4 + dummy1 + dummy2 + dummy3\nDesign Matrix:\n8×7 DataFrame\n│ Row │ factor1 │ factor2 │ factor3 │ factor4 │ dummy1 │ dummy2 │ dummy3 │\n│     │ Int64   │ Int64   │ Int64   │ Int64   │ Int64  │ Int64  │ Int64  │\n├─────┼─────────┼─────────┼─────────┼─────────┼────────┼────────┼────────┤\n│ 1   │ 1       │ 1       │ 1       │ 1       │ 1      │ 1      │ 1      │\n│ 2   │ -1      │ 1       │ -1      │ 1       │ 1      │ -1     │ -1     │\n│ 3   │ 1       │ -1      │ 1       │ 1       │ -1     │ -1     │ -1     │\n│ 4   │ -1      │ 1       │ 1       │ -1      │ -1     │ -1     │ 1      │\n│ 5   │ 1       │ 1       │ -1      │ -1      │ -1     │ 1      │ -1     │\n│ 6   │ 1       │ -1      │ -1      │ -1      │ 1      │ -1     │ 1      │\n│ 7   │ -1      │ -1      │ -1      │ 1       │ -1     │ 1      │ 1      │\n│ 8   │ -1      │ -1      │ 1       │ -1      │ 1      │ 1      │ -1     │\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.PlackettBurman-Tuple{StatsModels.FormulaTerm}","page":"Public","title":"ExperimentalDesign.PlackettBurman","text":"PlackettBurman(formula::StatsModels.FormulaTerm; symbol_encoding) -> PlackettBurman\n\n\njulia> PlackettBurman(@formula(y ~ x1 + x2 + x3 + x4))\nPlackettBurman\nDimension: (8, 7)\nFactors: (:x1, :x2, :x3, :x4)\nDummy Factors: (:dummy1, :dummy2, :dummy3)\nFormula: y ~ -1 + x1 + x2 + x3 + x4 + dummy1 + dummy2 + dummy3\nDesign Matrix:\n8×7 DataFrame\n│ Row │ x1    │ x2    │ x3    │ x4    │ dummy1 │ dummy2 │ dummy3 │\n│     │ Int64 │ Int64 │ Int64 │ Int64 │ Int64  │ Int64  │ Int64  │\n├─────┼───────┼───────┼───────┼───────┼────────┼────────┼────────┤\n│ 1   │ 1     │ 1     │ 1     │ 1     │ 1      │ 1      │ 1      │\n│ 2   │ -1    │ 1     │ -1    │ 1     │ 1      │ -1     │ -1     │\n│ 3   │ 1     │ -1    │ 1     │ 1     │ -1     │ -1     │ -1     │\n│ 4   │ -1    │ 1     │ 1     │ -1    │ -1     │ -1     │ 1      │\n│ 5   │ 1     │ 1     │ -1    │ -1    │ -1     │ 1      │ -1     │\n│ 6   │ 1     │ -1    │ -1    │ -1    │ 1      │ -1     │ 1      │\n│ 7   │ -1    │ -1    │ -1    │ 1     │ -1     │ 1      │ 1      │\n│ 8   │ -1    │ -1    │ 1     │ -1    │ 1      │ 1      │ -1     │\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#Base.rand","page":"Public","title":"Base.rand","text":"rand(distribution::CategoricalFactor) -> Any\nrand(distribution::CategoricalFactor, n::Int64) -> Any\n\n\njulia> rand(CategoricalFactor([:a, :b, 2, 1.0]), 6)\n6-element Array{Any,1}:\n  :b\n 2\n 1.0\n 2\n  :a\n 2\n\n\n\n\n\n\n","category":"function"},{"location":"lib/public/#Base.rand","page":"Public","title":"Base.rand","text":"rand(distribution::DesignDistribution) -> ExperimentalDesign.RandomDesign\nrand(distribution::DesignDistribution, n::Int64) -> ExperimentalDesign.RandomDesign\n\n\njulia> rand(DesignDistribution((f1 = Uniform(2, 3), f2 = DiscreteUniform(-1, 5), f3 = Uniform(5, 10))), 12)\nExperimentalDesign.RandomDesign\nDimension: (12, 3)\nFactors: (f1 = Distributions.Uniform{Float64}(a=2.0, b=3.0), f2 = Distributions.DiscreteUniform(a=-1, b=5), f3 = Distributions.Uniform{Float64}(a=5.0, b=10.0))\nFormula: 0 ~ f1 + f2 + f3\nDesign Matrix:\n12×3 DataFrame\n│ Row │ f1      │ f2      │ f3      │\n│     │ Float64 │ Float64 │ Float64 │\n├─────┼─────────┼─────────┼─────────┤\n│ 1   │ 2.04922 │ 5.0     │ 5.72029 │\n│ 2   │ 2.59117 │ 1.0     │ 8.48192 │\n│ 3   │ 2.77148 │ 4.0     │ 5.42517 │\n│ 4   │ 2.25659 │ 5.0     │ 6.75815 │\n│ 5   │ 2.64968 │ 5.0     │ 8.94993 │\n│ 6   │ 2.31523 │ 4.0     │ 5.89413 │\n│ 7   │ 2.86526 │ 3.0     │ 9.8807  │\n│ 8   │ 2.44753 │ -1.0    │ 9.24986 │\n│ 9   │ 2.08284 │ 3.0     │ 9.99911 │\n│ 10  │ 2.81201 │ 2.0     │ 7.28848 │\n│ 11  │ 2.64232 │ 4.0     │ 5.34507 │\n│ 12  │ 2.08189 │ 4.0     │ 8.303   │\n\n\n\n\n\n\n","category":"function"},{"location":"lib/public/#ExperimentalDesign.d_criterion-Tuple{Any}","page":"Public","title":"ExperimentalDesign.d_criterion","text":"d_criterion(model_matrix::Any; tolerance) -> Any\n\n\njulia> 1 + 1\n2\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.explicit_fullfactorial-Tuple{Base.Iterators.ProductIterator}","page":"Public","title":"ExperimentalDesign.explicit_fullfactorial","text":"explicit_fullfactorial(iterator::Base.Iterators.ProductIterator) -> Any\n\n\nReceives  a  Base.Iterators.ProductIterator  and  computes  an  explicit  full factorial design. The generated array is exponentially large.\n\njulia> explicit_fullfactorial(fullfactorial(([-1, 1], [:a, :b, :c])))\n6×2 Array{Any,2}:\n -1  :a\n  1  :a\n -1  :b\n  1  :b\n -1  :c\n  1  :c\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.explicit_fullfactorial-Tuple{Tuple}","page":"Public","title":"ExperimentalDesign.explicit_fullfactorial","text":"explicit_fullfactorial(factors::Tuple) -> Any\n\n\nReceives a tuple of arrays  representing categorical factor levels, and computes an explicit full factorial design. The generated array is exponentially large.\n\njulia> explicit_fullfactorial(([-1, 1], [:a, :b, :c], [1, 2]))\n12×3 Array{Any,2}:\n -1  :a  1\n  1  :a  1\n -1  :b  1\n  1  :b  1\n -1  :c  1\n  1  :c  1\n -1  :a  2\n  1  :a  2\n -1  :b  2\n  1  :b  2\n -1  :c  2\n  1  :c  2\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.fullfactorial-Tuple{Tuple}","page":"Public","title":"ExperimentalDesign.fullfactorial","text":"fullfactorial(factors::Tuple) -> Base.Iterators.ProductIterator\n\n\nReceives a tuple of arrays representing categorical factor levels, and returns a Base.Iterators.ProductIterator.   This allows  full  factorial  designs to  be arbitrarily large and  only be computed as needed.  To  compute an explicit full factorial design, use explicit_fullfactorial.\n\njulia> fullfactorial(Tuple([-1, 1] for i = 1:10))\nBase.Iterators.ProductIterator{NTuple{10,Array{Int64,1}}}(([-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1], [-1, 1]))\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.isplackettburman-Tuple{Array{Int64,2}}","page":"Public","title":"ExperimentalDesign.isplackettburman","text":"isplackettburman(d::Array{Int64,2}) -> Bool\n\n\nTo check if  a given design is  a Plackett-Burman design, we must check for the following properties, obtained in the original Plackett-Burman paper:\n\nEach component is replicated at each  of its values the same number of times, that is, the sum of elements in each column is zero\nEach  pair of components  occur together at  every combination of  values the same number of times, that is, the sum of each pair of columns will produce a column with the  same number of occurrences  of 2 and -2,  and twice that number of occurrences of 0\n\nPlackett,  R.L. and  Burman, J.P.,  1946. The  design of  optimum multifactorial experiments. Biometrika, 33(4), pp.305-325.\n\njulia> isplackettburman(plackettburman(2))\ntrue\n\njulia> isplackettburman(plackettburman(4))\ntrue\n\njulia> isplackettburman(plackettburman(16))\ntrue\n\njulia> isplackettburman(rand(Int, 4,4))\nfalse\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.kl_exchange-Tuple{StatsModels.FormulaTerm,DataFrames.DataFrame}","page":"Public","title":"ExperimentalDesign.kl_exchange","text":"kl_exchange(formula::StatsModels.FormulaTerm, candidate_set::DataFrames.DataFrame; tolerance, seed_design_size, max_iterations, experiments, design_k, candidates_l)\n\n\nBuilds optimum designs using the KL exchange algorithm, as described by Atkinson et al..   The ideia is  to iteratively swap  K design elements  with L candidate elements,  in order to  maximize an optimality  criterion.  Although any criteria based  on information matrices could be used  here, the KL exchange algorithm     leverages    determinant     properties     to    optimize     the d_criterion.\n\nAtkinson, A., Donev, A., & Tobias, R. (2007). Optimum experimental designs, with SAS (Vol. 34). Oxford University Press, Chapter 12.\n\njulia> candidates = FullFactorial(fill([-1, 0, 1], 5));\n\njulia> candidates_formula = ConstantTerm(0) ~ sum(Term.(Symbol.(names(candidates.matrix))));\n\njulia> candidates_formula = FormulaTerm(candidates_formula.lhs,\n                                        candidates_formula.rhs +\n                                        (@formula 0 ~ factor3 ^ 2).rhs);\n\njulia> selected_rows = kl_exchange(candidates_formula,\n                                   candidates.matrix,\n                                   seed_design_size = 2,\n                                   experiments = 11,\n                                   design_k = 11,\n                                   candidates_l = size(candidates.matrix, 1) - 11);\n\njulia> d_criterion(selected_rows)\n0.730476820204009\n\njulia> candidates.matrix[selected_rows.indices[1], :]\n11×5 DataFrames.DataFrame\n│ Row │ factor1 │ factor2 │ factor3 │ factor4 │ factor5 │\n│     │ Int64   │ Int64   │ Int64   │ Int64   │ Int64   │\n├─────┼─────────┼─────────┼─────────┼─────────┼─────────┤\n│ 1   │ -1      │ 1       │ -1      │ 1       │ 1       │\n│ 2   │ 1       │ -1      │ 1       │ 1       │ 1       │\n│ 3   │ 1       │ 1       │ 0       │ -1      │ 1       │\n│ 4   │ -1      │ 1       │ 1       │ -1      │ 1       │\n│ 5   │ -1      │ -1      │ 0       │ 1       │ 1       │\n│ 6   │ 1       │ 1       │ 1       │ 1       │ -1      │\n│ 7   │ -1      │ -1      │ 1       │ -1      │ -1      │\n│ 8   │ 1       │ -1      │ 0       │ 1       │ -1      │\n│ 9   │ -1      │ 1       │ -1      │ 1       │ -1      │\n│ 10  │ -1      │ 1       │ 0       │ -1      │ -1      │\n│ 11  │ 1       │ -1      │ -1      │ -1      │ -1      │\n\njulia> n_candidates = 3000;\n\njulia> n_factors = 30;\n\njulia> design_generator = DesignDistribution(Distributions.Uniform(0, 1), n_factors);\n\njulia> candidates = rand(design_generator, n_candidates);\n\njulia> selected_rows = kl_exchange(candidates.formula,\n                                   candidates.matrix,\n                                   experiments = n_factors + 2,\n                                   design_k = n_factors,\n                                   candidates_l = round(Int, (n_candidates - n_factors) / 2));\n\njulia> d_criterion(selected_rows)\n0.07904918814259465\n\njulia> candidates.matrix[selected_rows.indices[1], :]\n32×30 DataFrame. Omitted printing of 24 columns\n│ Row │ factor1    │ factor2   │ factor3   │ factor4   │ factor5   │ factor6   │\n│     │ Float64    │ Float64   │ Float64   │ Float64   │ Float64   │ Float64   │\n├─────┼────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤\n│ 1   │ 0.832487   │ 0.875035  │ 0.433788  │ 0.805738  │ 0.481517  │ 0.382364  │\n│ 2   │ 0.761498   │ 0.0398811 │ 0.998741  │ 0.632286  │ 0.29121   │ 0.55497   │\n│ 3   │ 0.639793   │ 0.49432   │ 0.678795  │ 0.685274  │ 0.269203  │ 0.896397  │\n│ 4   │ 0.00370418 │ 0.5139    │ 0.0630861 │ 0.175531  │ 0.570244  │ 0.60512   │\n│ 5   │ 0.698514   │ 0.207593  │ 0.766977  │ 0.0719112 │ 0.665055  │ 0.499193  │\n│ 6   │ 0.412744   │ 0.785263  │ 0.432861  │ 0.909976  │ 0.642545  │ 0.832013  │\n│ 7   │ 0.353518   │ 0.116685  │ 0.140352  │ 0.929095  │ 0.0484366 │ 0.838524  │\n⋮\n│ 25  │ 0.577343   │ 0.790529  │ 0.763994  │ 0.548131  │ 0.402066  │ 0.293336  │\n│ 26  │ 0.0941214  │ 0.52891   │ 0.293733  │ 0.92224   │ 0.982713  │ 0.929744  │\n│ 27  │ 0.0256128  │ 0.723072  │ 0.23927   │ 0.993221  │ 0.521335  │ 0.407319  │\n│ 28  │ 0.221631   │ 0.375189  │ 0.67251   │ 0.658613  │ 0.272448  │ 0.230501  │\n│ 29  │ 0.843885   │ 0.375497  │ 0.116069  │ 0.989937  │ 0.935444  │ 0.0466424 │\n│ 30  │ 0.181611   │ 0.956062  │ 0.98997   │ 0.80352   │ 0.994557  │ 0.203797  │\n│ 31  │ 0.738407   │ 0.25753   │ 0.929704  │ 0.407498  │ 0.0934688 │ 0.99359   │\n│ 32  │ 0.228975   │ 0.81272   │ 0.803107  │ 0.931826  │ 0.598371  │ 0.588823  │\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.next_offset_divisible_prime-NTuple{4,Int64}","page":"Public","title":"ExperimentalDesign.next_offset_divisible_prime","text":"next_offset_divisible_prime(n::Int64, offset::Int64, divisor::Int64, tries::Int64) -> Int64\n\n\nGets the next prime p, starting from n, for which (p + offset) % divisor == 0 holds.\n\njulia> next_offset_divisible_prime(3, 1, 4, 1000)\n3\n\njulia> next_offset_divisible_prime(5, 1, 4, 1000)\n7\n\njulia> next_offset_divisible_prime(4, 1, 4, 1000)\n7\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.paley-Tuple{Array{Int64,2}}","page":"Public","title":"ExperimentalDesign.paley","text":"paley(matrix::Array{Int64,2}) -> Array{Int64,2}\n\n\nThe Paley  construction is a method for constructing Hadamard matrices using finite fields.\n\njulia> paley(Matrix{Int}(undef, 8, 8))\n8×8 Array{Int64,2}:\n -1   1   1  -1   1   1   1  -1\n  1   1  -1   1   1   1  -1  -1\n  1  -1   1   1   1  -1  -1   1\n -1   1   1   1  -1  -1   1   1\n  1   1   1  -1  -1   1   1  -1\n  1   1  -1  -1   1   1  -1   1\n  1  -1  -1   1   1  -1   1   1\n -1  -1   1   1  -1   1   1   1\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.plackettburman-Tuple{Int64}","page":"Public","title":"ExperimentalDesign.plackettburman","text":"plackettburman(matrix_size::Int64) -> Array{Int64,2}\n\n\nConstructs a Plackett-Burman  design with size matrix_size if  possible, or to the closest, largest, number for which it is possible.\n\njulia> plackettburman(4)\n8×7 Array{Int64,2}:\n  1   1   1   1   1   1   1\n -1   1  -1   1   1  -1  -1\n  1  -1   1   1  -1  -1  -1\n -1   1   1  -1  -1  -1   1\n  1   1  -1  -1  -1   1  -1\n  1  -1  -1  -1   1  -1   1\n -1  -1  -1   1  -1   1   1\n -1  -1   1  -1   1   1  -1\n\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#ExperimentalDesign.random_design!-Tuple{Array{Float64,2},Tuple,Int64}","page":"Public","title":"ExperimentalDesign.random_design!","text":"random_design!(design::Array{Float64,2}, distributions::Tuple, n::Int64) -> Array{Float64,2}\n\n\njulia> random_design!(zeros(10, 3), (Uniform(2, 3), DiscreteUniform(-1, 5), Uniform(5, 10)), 10)\n10×3 Array{Float64,2}:\n 2.04922  -1.0  6.44508\n 2.59117   3.0  5.57183\n 2.77148   5.0  8.72759\n 2.25659   1.0  9.75865\n 2.64968   4.0  5.72029\n 2.31523   5.0  8.48192\n 2.86526   5.0  5.42517\n 2.44753   4.0  6.75815\n 2.08284   3.0  8.94993\n 2.81201  -1.0  5.89413\n\n\n\n\n\n\n","category":"method"},{"location":"#ExperimentalDesign.jl-1","page":"Home","title":"ExperimentalDesign.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Build Docs Test Coverage\n(Image: Build Status) (Image: ) (Image: Coverage Status) (Image: codecov.io)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"ExperimentalDesign.jl  provides  tools  for  Design  of  Experiments  in  Julia, enabling the construction  of designs for screening,  modeling, exploration, and optimization.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Development  on this  package  is  ongoing, so  expect  things  to change.  Pull requests are more than welcome!  Current features are:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Designs supporting categorical and continuous factors\nIntegration with StatsModels @formula\nFull factorial designs\nPlackett-Burman designs for screening\nFlexible random designs using the Distributions package\nVariance-optimizing designs for several criteria","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Intended features include the ones provided by R packages such as DoE.base, FrF2, and AlgDesign.","category":"page"},{"location":"#Library-Outline-1","page":"Home","title":"Library Outline","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\"lib/public.md\", \"lib/internals.md\"]","category":"page"},{"location":"#main-index-1","page":"Home","title":"Index","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\"lib/public.md\", \"lib/internals.md\"]","category":"page"}]
}