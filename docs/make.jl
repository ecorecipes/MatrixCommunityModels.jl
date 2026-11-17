using Documenter
using MatrixCommunityModels

# Tutorials are the rendered vignettes, copied here by scripts/sync_vignettes.jl.
function tutorial_pages()
    dir = joinpath(@__DIR__, "src", "tutorials")
    isdir(dir) || return Pair{String,String}[]
    files = sort(filter(f -> endswith(f, ".md"), readdir(dir)))
    return map(files) do f
        title = f
        for line in eachline(joinpath(dir, f))
            m = match(r"^#\s+(.*)", line)
            if m !== nothing
                title = String(strip(m.captures[1]))
                break
            end
        end
        return title => "tutorials/" * f
    end
end

# The open communities (`species`, `compose`, ...) are defined in CategoricalPopulationDynamics.jl
# and re-exported; their docstrings are listed in api/communities.md. That module is passed to
# `modules` so that `@docs` finds them, and ignored by `checkdocs`, which covers this package's
# own exports. Their source links need the remote of that package when it is installed from
# GitHub rather than checked out (its commit is then unknown, so links point at `main`).
const CPD = MatrixCommunityModels.CategoricalPopulationDynamics
const CPD_DIR = pkgdir(CPD)
const CPD_REMOTE = Documenter.Remotes.GitHub("ecorecipes", "CategoricalPopulationDynamics.jl")
const REMOTES = isdir(joinpath(CPD_DIR, ".git")) || isfile(joinpath(CPD_DIR, ".git")) ? Dict() :
                Dict(CPD_DIR => (CPD_REMOTE, "main"))

makedocs(;
         modules = [MatrixCommunityModels, CPD],
         checkdocs_ignored_modules = [CPD],
         remotes = REMOTES,
         authors = "Simon Frost",
         sitename = "MatrixCommunityModels.jl",
         checkdocs = :exports,
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://ecorecipes.github.io/MatrixCommunityModels.jl",
                                  edit_link = "main",
                                  size_threshold = 300 * 1024),
         pages = ["Home" => "index.md",
                  "Tutorials" => tutorial_pages(),
                  "API Reference" => ["Composing communities" => "api/communities.md",
                                      "Analysis" => "api/analysis.md",
                                      "IPMs and kernels" => "api/kernels.md",
                                      "Published models" => ["Rivers and orchids" => "api/models_rivers.md",
                                                             "Bees and flowers" => "api/models_bees.md",
                                                             "Plants and trees" => "api/models_plants.md",
                                                             "Reefs and seas" => "api/models_marine.md",
                                                             "Competitors, predators and pathogens" => "api/models_animals.md"]]])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/MatrixCommunityModels.jl.git", push_preview = true)
end
