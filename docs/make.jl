using Documenter
using RCDesignSuite

makedocs(
    sitename    = "RCDesignSuite.jl",
    format      = Documenter.HTML(),
    modules     = [RCDesignSuite],
    repo        = "https://github.com/BYU-Aeronautics-Club/RCDesignSuite.jl/blob/{commit}{path}#L{line}",
    authors     = "Judd Mehr <juddmehr@gmail.com>",
)

deploydocs(
    repo = "github.com/BYU-Aeronautics-Club/RCDesignSuite.jl.git"
)
