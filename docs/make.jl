using HyperplanesAndMatroids
using Documenter

DocMeta.setdocmeta!(HyperplanesAndMatroids, :DocTestSetup, :(using HyperplanesAndMatroids); recursive=true)

makedocs(;
    modules=[HyperplanesAndMatroids],
    authors="lienkaemper <caitlienk@gmail.com> and contributors",
    repo="https://github.com/lienkaemper/HyperplanesAndMatroids.jl/blob/{commit}{path}#{line}",
    sitename="HyperplanesAndMatroids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lienkaemper.github.io/HyperplanesAndMatroids.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lienkaemper/HyperplanesAndMatroids.jl",
    devbranch="main",
)
