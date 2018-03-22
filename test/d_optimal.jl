using ExperimentalDesign

ranges = [3:3, 4:4, 5:5, 6:6, 7:7, 8:8, 9:9, 10:10, 11:11, 12:12, 5:16]
factors = [[Array{Float64, 1}(1:2) for i = 1:4] for k in 1:length(ranges)]
designs = 8000

scaling!(design, factors) = scale_boxdraper_encoding!(design, factors, scale_denominator = false)

sampled_subsets = sample_subsets(factors, ranges, designs, scale = scaling!)
