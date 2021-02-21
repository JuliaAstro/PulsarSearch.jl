using PulsarSearch
using Documenter

makedocs(;
    modules=[PulsarSearch],
    authors="Matteo Bachetti <matteo@matteobachetti.it> and contributors",
    repo="https://github.com/JuliaAstro/PulsarSearch.jl/blob/{commit}{path}#L{line}",
    sitename="PulsarSearch.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/PulsarSearch.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAstro/PulsarSearch.jl",
    devbranch = "main"
)
