# Slow checks, left out of the main suite for time: long runs, ensembles and the comparisons that need them. They
# run after the main suite when the environment variable MCM_LONG_TESTS is set to true:
#
#     MCM_LONG_TESTS=true julia --project -e 'using Pkg; Pkg.test()'
#
# and every week on GitHub Actions (.github/workflows/LongTests.yml). Tolerances follow the rule of the main suite
# (about 100 times the achieved agreement where only rounding limits it, 1.5 to 3 times where the reference program's
# own numerics do, printed precision for printed or published values), loosened to about 1e-5 relative where a long
# nonlinear run could differ between platforms by more than its achieved agreement.
using Test
using MatrixCommunityModels

# a reference file of test/reference, as rows of fields, without its header
long_reference(file) = [split(l, ",") for l in readlines(joinpath(@__DIR__, "..", "reference", file))[2:end]]
# the statistics the vignettes print, without the Statistics package (not a test dependency)
long_mean(x) = sum(x) / length(x)
function long_quantile(x, p)            # R's and Statistics.quantile's default (type 7)
    s = sort(collect(float.(x)))
    h = (length(s) - 1) * p + 1
    lo = floor(Int, h)
    return lo == length(s) ? s[end] : s[lo] + (h - lo) * (s[lo + 1] - s[lo])
end
long_median(x) = long_quantile(x, 0.5)

@testset "MatrixCommunityModels, long tests" begin
    include("test_grasslands_long.jl")
    include("test_forests_long.jl")
    include("test_sagebrush_long.jl")
    include("test_chytrid_long.jl")
    include("test_thistles_long.jl")
end
