using ExperimentalDesign, DataStructures

ranges = [8:16]
factors = [OrderedDict([(:f1, [1., 2.]),
                        (:f2, [1., 2.]),
                        (:f3, [1., 2.]),
                        (:f4, [1., 2.]),
                        (:f5, [1., 2.]),
                        (:f6, [1., 2.]),
                        (:f7, [1., 2.])
                       ]) for i in ranges]
designs = 8000

sampled_subsets = sample_subsets(factors, ranges, designs, scale = scale_boxdraper_encoding!)
