using MRIBuilder
using Documenter

DocMeta.setdocmeta!(MRIBuilder, :DocTestSetup, :(using MRIBuilder); recursive=true)

makedocs(;
    modules=[MRIBuilder],
    authors="Michiel Cottaar <Michiel.Cottaar@ndcn.ox.ac.uk>",
    sitename="MRIBuilder.jl",
    format=Documenter.HTML(;
        canonical="https://MichielCottaar.gitlab.io/MRIBuilder.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
