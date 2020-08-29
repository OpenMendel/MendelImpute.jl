using Documenter
using MendelImpute

makedocs(
    sitename = "MendelImpute",
    format = Documenter.HTML(),
    modules = [MendelImpute]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
