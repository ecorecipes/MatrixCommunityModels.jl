using Test
using MatrixCommunityModels

@testset "MatrixCommunityModels" begin
    include("test_verde_fish.jl")
    include("test_adapt3.jl")
    include("test_riparian.jl")
    include("test_analysis.jl")
    include("test_bees.jl")
    include("test_docstrings.jl")
end
