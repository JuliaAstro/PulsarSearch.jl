using psrsearch
using Documenter

makedocs(;
    modules=[psrsearch],
    authors="Matteo Bachetti <matteo@matteobachetti.it> and contributors",
    repo="https://github.com/JuliaAstro/psrsearch.jl/blob/{commit}{path}#L{line}",
    sitename="psrsearch.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/psrsearch.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAstro/psrsearch.jl",
    devbranch = "main"
)
