using ExperimentalDesign, DataStructures

ranges = [k:k for k = 20:1:24]

factors = OrderedDict{Symbol, Any}()

for i = 1:10
    factors[Symbol(:f, i)] = [:a, :b, :c]
end

factor_list = [factors for i in ranges]

designs = 4000

sampled_subsets = sample_subsets(factor_list, ranges, designs, check_bounds = true)
