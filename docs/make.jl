using Documenter
using RCDesignSuite

makedocs(
    modules     = [RCDesignSuite],
    sitename    = "RCDesignSuite.jl",
    format      = Documenter.HTML(),
    pages       = [
        "Intro" => "index.md",
        "Quick Start" => "tutorial.md",
        "Guided Examples" => "howto.md",
        "API Reference" => "reference.md",
        "Theory" => "theory.md"
    ],
    repo        = "https://github.com/BYU-Aeronautics-Club/RCDesignSuite.jl/blob/{commit}{path}#L{line}",
    authors     = "Judd Mehr <juddmehr@gmail.com>",
)

deploydocs(
    repo = "github.com/BYU-Aeronautics-Club/RCDesignSuite.jl.git"
)
