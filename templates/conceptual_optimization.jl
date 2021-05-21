#=

A first attempt at a conceptual  level optimization.

Authors: Judd Mehr,

=#

#############################################
##########          SETUP          ##########
#############################################
using SNOW

include("templates/previous_year_objectives.jl")

x0, lb, ub, p, c = setup2021()

function f!(con, x)

    J = obj2021(x,p,c)

    con = con2021(x,p,c)

    return J

end

nc = 4
lc = -Inf*ones(nc)  # lower bounds on constraints
uc = zeros(nc)  # upper bounds on constraints
options = Options(solver=IPOPT())  # choosing IPOPT solver

xopt, fopt, info = minimize(f!, x0, nc, lb, ub, lc, uc, options)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)