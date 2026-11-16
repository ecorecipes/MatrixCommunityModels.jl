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
    include("test_forests.jl")
    include("test_tabulated.jl")
    include("test_seafans.jl")
    include("test_sagebrush.jl")
    include("test_lemans.jl")
    include("test_solvers.jl")
    include("test_docstrings.jl")
end
