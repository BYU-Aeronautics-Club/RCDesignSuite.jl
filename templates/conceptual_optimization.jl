#=

A first attempt at a conceptual  level optimization.

Authors: Judd Mehr,

=#

#############################################
##########          SETUP          ##########
#############################################
using SNOW
using RCDesignSuite

include("previous_year_objectives.jl")

x0, lx, ux, p, c = setup2021()

function f!(con, x)

    J = obj2021(x,p,c)

    con2021!(con, x, p, c)

    return J

end

ng = 7
lg = -Inf*ones(ng)  # lower bounds on constraints
ug = zeros(ng)  # upper bounds on constraints
options = Options(solver=IPOPT())  # choosing IPOPT solver

xopt, fopt, info = minimize(f!, copy(x0), ng, lx, ux, lg, ug, options)

con = zeros(ng)
con2021!(con, xopt,p,c)
display(con)

println("xstar = ")
display(xopt)
println("fstar = ")
display(fopt)
println("info = ")
display(info)