using Documenter
using MendelImpute

makedocs(
    sitename = "MendelImpute",
    format = Documenter.HTML(),
    modules = [MendelImpute],
    pages = [
        "Home" => "index.md",
        "Phasing and Imputation" => "man/Phasing_and_Imputation.md",
        "Performance Gotchas" => "man/performance.md",
        "Estimating ancestry" => "man/painting.md",
        "Ultra compression" => "man/ultra+compress.md",
        "Run as script" => "man/script.md",
        "API" => "man/api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/OpenMendel/MendelImpute.jl.git",
    target = "build"
)
