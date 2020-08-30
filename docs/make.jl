using Documenter
using MendelImpute

makedocs(
    sitename = "MendelImpute",
    format = Documenter.HTML(),
    modules = [MendelImpute]
    pages = [
        "Home" => "index.md",
        "API"  => "man/api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/biona001/MendelImpute.jl.git",
    target = "build"
)
