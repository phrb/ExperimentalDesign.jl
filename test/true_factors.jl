using ExperimentalDesign, DataStructures

factors = 10

#ranges = [k:k for k = ((2 * factors) + 1):1:((2 * factors) + 5)]
ranges = [k:k for k = 16:1:31]

factors = OrderedDict{Symbol, Any}()

for i = 1:5
    factors[Symbol(:f, i)] = Int[-1, 0, 1]
end

for i = 6:10
    factors[Symbol(:f, i)] = [-1., 1.]
end

factor_list = [factors for i in ranges]

designs = 8000

sampled_subsets = sample_subsets(factor_list, ranges, designs, check_bounds = true)
