using MRIBuilder
using Documenter

DocMeta.setdocmeta!(MRIBuilder, :DocTestSetup, :(using MRIBuilder); recursive=true)

makedocs(;
    modules=[MRIBuilder],
    authors="Michiel Cottaar <Michiel.Cottaar@ndcn.ox.ac.uk>",
    repo=Remotes.GitLab("git.fmrib.ox.ac.uk", "ndcn0236", "MRIBuilder.jl"),
    sitename="MRIBuilder.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        size_threshold_ignore=["api.md"],
        description="Documentation for MRIBuilder.jl: an MRI sequence builder in Julia",
        footer=nothing,
        canonical="https://open.win.ox.ac.uk/pages/ndcn0236/mribuilder.jl/stable/",
    ),
    pages=[
        "Home" => "index.md",
        "Optimisation" => "sequence_optimisation.md",
        "Implemented sequences" => "implemented_sequences.md",
        "Scanners" => "scanners.md",
        "Internal API" => "api.md",
    ],
    warnonly=Documenter.except(:example_block),
)

if get(ENV, "CI_COMMIT_REF_NAME", "") == "main" || length(get(ENV, "CI_COMMIT_TAG", "")) > 0
    deploydocs(repo="git.fmrib.ox.ac.uk:ndcn0236/mribuilder.jl.git", branch="pages", devbranch="main")
else
    println("Skipping deployment, because we are local or on a secondary branch.")
end