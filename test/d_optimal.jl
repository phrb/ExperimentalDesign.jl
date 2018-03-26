using ExperimentalDesign, DataStructures

factors = 10

ranges = [k:k for k = factors:2:(2 * factors)]

factor_dict = [OrderedDict([(Symbol(:f, j), [1., 2.]) for j = 1:factors]) for i in ranges]

designs = 400

sampled_subsets = sample_subsets(factor_dict, ranges, designs, check_bounds = false)
