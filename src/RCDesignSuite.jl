module RCDesignSuite

#import packages
import PyPlot; plt = PyPlot
import DataFrames
import CSV
import Xfoil
using Printf

#include code
include("aircrafttypes.jl")
include("commonplots.jl")
include("conceptual.jl")
include("../templates/sensitivity.jl")

end # module
