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

makedocs(;
         modules = [MatrixCommunityModels],
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
                  "API Reference" => "api.md"])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/MatrixCommunityModels.jl.git", push_preview = true)
end
