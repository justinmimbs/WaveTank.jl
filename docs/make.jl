using Documenter
using WaveTank

DocMeta.setdocmeta!(WaveTank, :DocTestSetup, :(using WaveTank); recursive=true)

makedocs(
    sitename="WaveTank.jl",
    modules=[WaveTank],
    format=Documenter.HTML(prettyurls=false),
)
