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
    include("test_annualplants.jl")
    include("test_waterfowl.jl")
    include("test_solvers.jl")
    include("test_resilience.jl")
    include("test_management.jl")
    include("test_docstrings.jl")
end

# the slow checks (test/long), only when asked: MCM_LONG_TESTS=true
if lowercase(get(ENV, "MCM_LONG_TESTS", "false")) in ("true", "1", "yes")
    include(joinpath("long", "runtests.jl"))
end

# The DecisionAnalysis.jl adapter (ext/MatrixCommunityModelsDecisionAnalysisExt.jl) is tested in an
# environment of its own, test/decision_analysis, in a second process: DecisionAnalysis.jl needs
# Julia 1.12, so it cannot be a test dependency of a package that supports 1.11. The tests are
# skipped below 1.12, with MCM_DECISION_TESTS=false, and (as a broken test, with a warning) when
# DecisionAnalysis.jl cannot be installed.
if VERSION >= v"1.12" && lowercase(get(ENV, "MCM_DECISION_TESTS", "true")) in ("true", "1", "yes")
    @testset "DecisionAnalysis adapter (own environment)" begin
        project = joinpath(@__DIR__, "decision_analysis")
        julia = Base.julia_cmd()
        env = filter(p -> !startswith(first(p), "JULIA_LOAD_PATH") && first(p) != "JULIA_PROJECT", collect(ENV))
        run_in(cmd) = success(pipeline(setenv(cmd, env); stdout, stderr))
        if run_in(`$julia --project=$project -e "using Pkg; Pkg.instantiate()"`)
            @test run_in(`$julia --project=$project $(joinpath(project, "runtests.jl"))`)
        else
            @warn "DecisionAnalysis.jl could not be installed; its adapter was not tested"
            @test_skip false
        end
    end
else
    @info "The DecisionAnalysis.jl adapter was not tested" VERSION get(ENV, "MCM_DECISION_TESTS", "true")
end
