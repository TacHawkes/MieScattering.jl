using Documenter
using MieScattering

makedocs(
    modules = [MieScattering],
    authors = "Oliver Kliebisch <oliver@kliebisch.net> and contributors",
    sitename = "MieScattering.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        assets = String[],
    ),
    pages = [
        "Introduction" => "index.md",
        "Performance" => "performance.md",
        "API" => "api.md",
    ],
)

deploydocs(repo = "github.com/TacHawkes/MieScattering.jl.git")
