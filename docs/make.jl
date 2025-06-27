using Documenter
using GalaxyProfiles
# using GalaxyProfiles.SurfaceDensities

# The `format` below makes it so that urls are set to "pretty" if you are pushing them to a hosting service, and basic if you are just using them locally to make browsing easier.

DocMeta.setdocmeta!(GalaxyProfiles, :DocTestSetup, :(using GalaxyProfiles); recursive=true)

makedocs(
    sitename = "GalaxyProfiles.jl",
    modules = [GalaxyProfiles],
    format = Documenter.HTML(;prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Chris Garling",
    pages = ["index.md","types.md","methods.md","guide.md","units.md","docindex.md"],
    push_preview = true,
    doctest = true,
    linkcheck = true,
    # Do not error if we are missing a docstring in the module or if an external link is invalid
    warnonly = [:missing_docs, :linkcheck]    
)

# deploydocs(;
#     repo = "github.com/cgarling/GalaxyProfiles.jl.git",
#     versions = ["stable" => "v^", "v#.#"],
#     push_preview=true,
# )
