using Documenter
using StarStuff

makedocs(
    sitename = "StarStuff",
    format = Documenter.HTML(),
    modules = [StarStuff],
    pages = [
        "Random Number Generator" => "generator.md",
        "Particle Trees" => "particletree.md"
    ]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
