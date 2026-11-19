using Test
using MatrixCommunityModels

@testset "MatrixCommunityModels" begin
    include("test_verde_fish.jl")
    include("test_adapt3.jl")
    include("test_riparian.jl")
    include("test_analysis.jl")
    include("test_kernels.jl")
    include("test_bees.jl")
    include("test_thistles.jl")
    include("test_recruitment.jl")
    include("test_forests.jl")
    include("test_tabulated.jl")
    include("test_seafans.jl")
    include("test_sagebrush.jl")
    include("test_lemans.jl")
    include("test_corals.jl")
    include("test_grasslands.jl")
    include("test_tribolium.jl")
    include("test_flowregimes.jl")
    include("test_predatorprey.jl")
    include("test_chytrid.jl")
    include("test_mizer.jl")
    include("test_solvers.jl")
    include("test_resilience.jl")
    include("test_docstrings.jl")
end

# the slow checks (test/long), only when asked: MCM_LONG_TESTS=true
if lowercase(get(ENV, "MCM_LONG_TESTS", "false")) in ("true", "1", "yes")
    include(joinpath("long", "runtests.jl"))
end
