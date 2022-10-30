using Documenter
using MieScattering

makedocs(modules = [MieScattering],
         authors = "Oliver Kliebisch <oliver@kliebisch.net> and contributors",
         repo = "https://github.com/tachawkes/miescattering.jl/blob/{commit}{path}#L{line}",
         sitename = "MieScattering.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  assets = String[]),
         pages = [
             "Introduction => "index.md",
         ])

deploydocs(repo = "github.com/TacHawkes/MieScattering.jl.git")
