#=

A first attempt at a conceptual  level optimization.

Authors: Judd Mehr,

=#

#############################################
##########          SETUP          ##########
#############################################
using SNOW
using RCDesignSuite

include("templates/previous_year_objectives.jl")

x0, lx, ux, p, c = setup2021()

function f!(con, x)

    J = obj2021(x,p,c)

    con .= con2021(x,p,c)

    return J

end

ng = 4
lg = -Inf*ones(ng)  # lower bounds on constraints
ug = zeros(ng)  # upper bounds on constraints
options = Options(solver=IPOPT())  # choosing IPOPT solver

xopt, fopt, info = minimize(f!, copy(x0), ng, lx, ux, lg, ug, options)

con2021(xopt,p,c)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)