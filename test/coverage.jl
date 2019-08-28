# Only run coverage from linux build on travis.
get(ENV, "TRAVIS_OS_NAME", "")       == "linux"   || exit()
using Pkg
Pkg.add("Coverage")
using Coverage

cd(joinpath(dirname(@__FILE__), "..")) do
    Codecov.submit(Codecov.process_folder())
    Coveralls.submit(Coveralls.process_folder())
end
