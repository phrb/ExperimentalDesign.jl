var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ExperimentalDesign.d_efficiency_lower_bound-Tuple{Array{Float64,2}}",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.d_efficiency_lower_bound",
    "category": "method",
    "text": "d_efficiency_lower_bound(model_matrix::Array{Float64, 2})\n\nCompute a lower bound for the D-efficiency of a given design\'s model matrix according to Castillo\'s \"Process Optimization : A Statistical Approach\".\n\njulia> using ExperimentalDesign, DataStructures, StatsModels\n\njulia> A = plackett_burman(4)\n4×3 Array{Int64,2}:\n  1   1   1\n  1  -1  -1\n -1  -1   1\n -1   1  -1\n\njulia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])\nDataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:\n  :f1 => [-1.0, 1.0]\n  :f2 => [-1.0, 1.0]\n  :f3 => [-1.0, 1.0]\n\njulia> m = generate_model_matrix(@formula(y ~ f1 + f2 + f3), float(A), factors)\n4×4 Array{Float64,2}:\n 1.0   1.0   1.0   1.0\n 1.0   1.0  -1.0  -1.0\n 1.0  -1.0  -1.0   1.0\n 1.0  -1.0   1.0  -1.0\n\njulia> d_efficiency_lower_bound(m)\n1.0\n\n\nFormula\n\nFor a design A_np with n experiments or rows, p factors or columns, and model matrix mathbfX, the lower bound for the D-efficiency of A is:\n\nD_eff^(L) = dfracmathbfX^primemathbfX^1pn\n\n\n\n"
},

{
    "location": "index.html#ExperimentalDesign.d_optimality-Tuple{Array{Float64,2}}",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.d_optimality",
    "category": "method",
    "text": "d_optimality(model_matrix::Array{Float64, 2})\n\nCompute the D-optimality of a given design\'s model matrix.\n\nExamples\n\njulia> using ExperimentalDesign, DataStructures, StatsModels\n\njulia> A = plackett_burman(4)\n4×3 Array{Int64,2}:\n  1   1   1\n  1  -1  -1\n -1  -1   1\n -1   1  -1\n\njulia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])\nDataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:\n  :f1 => [-1.0, 1.0]\n  :f2 => [-1.0, 1.0]\n  :f3 => [-1.0, 1.0]\n\njulia> m = generate_model_matrix(@formula(y ~ f1 + f2 + f3), float(A), factors)\n4×4 Array{Float64,2}:\n 1.0   1.0   1.0   1.0\n 1.0   1.0  -1.0  -1.0\n 1.0  -1.0  -1.0   1.0\n 1.0  -1.0   1.0  -1.0\n\njulia> d_optimality(m)\n256.0\n\n\nFormula\n\nThe D-optimality of a design is the determinant of the information matrix mathbfX^primemathbfX, where mathbfX is the model matrix of a design.\n\n\n\n"
},

{
    "location": "index.html#ExperimentalDesign.expand_factors-Tuple{DataStructures.OrderedDict}",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.expand_factors",
    "category": "method",
    "text": "expand_factors(factors::OrderedDict)\n\nReturn an OrderedDict with true factors expanded to [0., 1.] 2-level numerical factors.\n\nFor a true factor with n levels, creates n - 1 2-level factors. Only one of these new factors can be at level 1. at each row of a design, encoding each level of the original factor. The first level of the original factor is encoded by all new factors being at level 0..\n\nExamples\n\njulia> using ExperimentalDesign, DataStructures, StatsModels\n\njulia> A = OrderedDict(:A => [1., 2.], :B => [1, 2, 3, 4], :C => [\"A\", \"B\"], :D => [1.2, 2.2, 3.4], :E => [:a, :b, :c, :d])\nDataStructures.OrderedDict{Symbol,Any} with 5 entries:\n  :A => [1.0, 2.0]\n  :B => [1, 2, 3, 4]\n  :C => String[\"A\", \"B\"]\n  :D => [1.2, 2.2, 3.4]\n  :E => Symbol[:a, :b, :c, :d]\n\njulia> expand_factors(A)\nDataStructures.OrderedDict{Symbol,Any} with 7 entries:\n  :A   => [1.0, 2.0]\n  :B   => [1, 2, 3, 4]\n  :C_1 => [0.0, 1.0]\n  :D   => [1.2, 2.2, 3.4]\n  :E_1 => [0.0, 1.0]\n  :E_2 => [0.0, 1.0]\n  :E_3 => [0.0, 1.0]\n\n\n\n\n"
},

{
    "location": "index.html#ExperimentalDesign.generate_model_matrix-Tuple{StatsModels.Formula,DataFrames.DataFrame,Dict}",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.generate_model_matrix",
    "category": "method",
    "text": "generate_model_matrix(formula::Formula,\n                      design::Array{Float64, 2},\n                      factors::OrderedDict;\n                      scale::Function = scale_boxdraper_encoding!)\n\nGenerate a DataFrame with a scaled model matrix for a given formula, design and factors.\n\nAssumes that formula is a linear relationship between all the factors in factors.\n\nExamples\n\njulia> using ExperimentalDesign, DataStructures, StatsModels\n\njulia> A = plackett_burman(4)\n4×3 Array{Int64,2}:\n  1   1   1\n  1  -1  -1\n -1  -1   1\n -1   1  -1\n\njulia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])\nDataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:\n  :f1 => [-1.0, 1.0]\n  :f2 => [-1.0, 1.0]\n  :f3 => [-1.0, 1.0]\n\njulia> m = generate_model_matrix(@formula(y ~ f1 + f2 + f3), float(A), factors)\n4×4 Array{Float64,2}:\n 1.0   1.0   1.0   1.0\n 1.0   1.0  -1.0  -1.0\n 1.0  -1.0  -1.0   1.0\n 1.0  -1.0   1.0  -1.0\n\n\n\n\n"
},

{
    "location": "index.html#ExperimentalDesign.get_expanded_values-Tuple{DataStructures.OrderedDict}",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.get_expanded_values",
    "category": "method",
    "text": "get_expanded_values(factors::OrderedDict)\n\nReturn a mapping from an original level of a factor to a column in an expanded design. Does not expand factors.\n\nExamples\n\njulia> using ExperimentalDesign, DataStructures, StatsModels\n\njulia> A = OrderedDict(:A => [1., 2.], :B => [1, 2, 3, 4], :C => [\"A\", \"B\"], :D => [1.2, 2.2, 3.4], :E => [:a, :b, :c, :d])\nDataStructures.OrderedDict{Symbol,Any} with 5 entries:\n  :A => [1.0, 2.0]\n  :B => [1, 2, 3, 4]\n  :C => String[\"A\", \"B\"]\n  :D => [1.2, 2.2, 3.4]\n  :E => Symbol[:a, :b, :c, :d]\n\njulia> get_expanded_values(A)\nDataStructures.OrderedDict{Any,Any} with 2 entries:\n  :C => DataStructures.OrderedDict{Any,Symbol}(\"B\"=>:C_1)\n  :E => DataStructures.OrderedDict{Any,Symbol}(:b=>:E_1,:c=>:E_2,:d=>:E_3)\n\n\n\n\n"
},

{
    "location": "index.html#ExperimentalDesign.scale_boxdraper_encoding!-Union{Tuple{Array{Float64,2},Array{T,1}}, Tuple{T}} where T",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.scale_boxdraper_encoding!",
    "category": "method",
    "text": "scale_boxdraper_encoding!(design::Array{Float64, 2},\n                          factors::Array{T, 1};\n                          scale_denominator = true) where T <: Any\n\nScale factors of a design using the Box and Draper\'s coding convention from \"Response Surfaces, Mixtures, and Ridge Analyses\".\n\nExamples\n\nWith denominator scaling:\n\njulia> using ExperimentalDesign, DataStructures, StatsModels\n\njulia> A = float(plackett_burman(4))\n4×3 Array{Float64,2}:\n  1.0   1.0   1.0\n  1.0  -1.0  -1.0\n -1.0  -1.0   1.0\n -1.0   1.0  -1.0\n\njulia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])\nDataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:\n  :f1 => [-1.0, 1.0]\n  :f2 => [-1.0, 1.0]\n  :f3 => [-1.0, 1.0]\n\njulia> scale_boxdraper_encoding!(A, collect(values(factors)))\n4×3 Array{Float64,2}:\n  1.0   1.0   1.0\n  1.0  -1.0  -1.0\n -1.0  -1.0   1.0\n -1.0   1.0  -1.0\n\njulia> all([isapprox(1.0, sqrt(sum(A[:, i] .^ 2.0) / size(A, 1))) for i = 1:size(A, 2)])\ntrue\n\n\nWithout denominator scaling:\n\njulia> using ExperimentalDesign, DataStructures, StatsModels\n\njulia> A = float(plackett_burman(4))\n4×3 Array{Float64,2}:\n  1.0   1.0   1.0\n  1.0  -1.0  -1.0\n -1.0  -1.0   1.0\n -1.0   1.0  -1.0\n\njulia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])\nDataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:\n  :f1 => [-1.0, 1.0]\n  :f2 => [-1.0, 1.0]\n  :f3 => [-1.0, 1.0]\n\njulia> scale_boxdraper_encoding!(A, collect(values(factors)), scale_denominator = false)\n4×3 Array{Float64,2}:\n  0.5   0.5   0.5\n  0.5  -0.5  -0.5\n -0.5  -0.5   0.5\n -0.5   0.5  -0.5\n\njulia> all([isapprox(1.0, sqrt(sum(A[:, i] .^ 2.0))) for i = 1:size(A, 2)])\ntrue\n\n\nFormula\n\nFor a design D_np with n experiments or rows and p factors or columns, scale each factor mathbfx_i to mathbfx_i^s according to:\n\nmathbfx_i^s = dfracmathbfx_i - barmathbfx_iS_i\n\nWhere barmathbfx_i is mean of factor definition values in the factors parameter and:\n\nS_i^2 = dfrac1n sumlimits_j = 1^n(x_ij - barmathbfx_i)^2\n\nIf we pass scale_denominator = false, S_i becomes:\n\nS_i^2 = sumlimits_j = 1^n(x_ij - barmathbfx_i)^2\n\n\n\n"
},

{
    "location": "index.html#ExperimentalDesign.scale_orthogonal!-Union{Tuple{Array{Float64,2},Array{T,1}}, Tuple{T}} where T",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.scale_orthogonal!",
    "category": "method",
    "text": "scale_orthogonal!(design::Array{Float64, 2}, factors::Array{T, 1}) where T <: Any\n\nOrthogonally scale and center factors of a design using design and factor limits.\n\nExamples\n\njulia> using ExperimentalDesign, DataStructures, StatsModels\n\njulia> A = [5. 2. -1.; 5. 3. 0.; -5. 1. -2.]\n3×3 Array{Float64,2}:\n  5.0  2.0  -1.0\n  5.0  3.0   0.0\n -5.0  1.0  -2.0\n\njulia> factors = OrderedDict([(:f1, [-5., 0., 5.]), (:f2, [1., 2., 3.]), (:f3, [-2., -1., 0.])])\nDataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:\n  :f1 => [-5.0, 0.0, 5.0]\n  :f2 => [1.0, 2.0, 3.0]\n  :f3 => [-2.0, -1.0, 0.0]\n\njulia> scale_orthogonal!(A, collect(values(factors)))\n3×3 Array{Float64,2}:\n  1.0   0.0   0.0\n  1.0   1.0   1.0\n -1.0  -1.0  -1.0\n\njulia> A\n3×3 Array{Float64,2}:\n  1.0   0.0   0.0\n  1.0   1.0   1.0\n -1.0  -1.0  -1.0\n\n\nFormula\n\nFor a design D_np with n experiments or rows and p factors or columns, scale each factor mathbfx_i to mathbfx_i^s according to:\n\nmathbfx_i^s = dfracmathbfx_i - barMbarM_def\n\nWhere mathbfx_def is the factor defined in the factors parameter and:\n\nbarM = (max(mathbfx_i) + min(mathbfx_i))  2\n\nbarM_def = (max(mathbfx_def) - min(mathbfx_def))  2\n\n\n\n"
},

{
    "location": "index.html#ExperimentalDesign.jl-Documentation-1",
    "page": "ExperimentalDesign.jl Documentation",
    "title": "ExperimentalDesign.jl Documentation",
    "category": "section",
    "text": "Modules = [ExperimentalDesign]"
},

]}
