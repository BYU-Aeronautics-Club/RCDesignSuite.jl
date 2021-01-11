#=
Author: Judd Mehr

Several ready-made plotting functions for common plots
=#

export plotaf



#######################
####### General #######
#######################
"""
"""
function plotbasic()

end


"""
"""
function ploterrbars()

end


#######################
###### AIRFOILS #######
#######################
"""
    plotaf(x, y; width=7)

Plot airfoil(s) with unit chord.
"""
function plotaf(x, y, labels;
            width       = 7,
            colorcycle  = ["C0"; "C1"; "C3"; "C0"; "C1"; "C3"], stylecycle  = ["-";"-";"-";"--";"--";"--"]
            )

    if length(x[1,:]) != length(colorcycle)
        @error("You are trying to plot more than 6 airfoils, you need to define longer color and style cycles.")
    end

    height = width*2.1*maximum(y)
    figure(figsize=(width,height))
    for i=1:length(x[1,:])
        plot(x[i,:]./maximum(x[i,:]),y[i,:]./maximum(x[i,:]),color=colorcycle[i],style=stylecycle[i],label=labels[i])
    end

    return false

end


"""
"""
function plotclcd()

end


"""
"""
function plotliftdist()

end


#######################
###### STABILITY ######
#######################
"""
"""
function plotcm()

end


"""
"""
function plotrootlocus()

end


#######################
##### PERFORMANCE #####
#######################
"""
"""
function plotpropefficiency()

end


"""
"""
function plotdragbars()

end


"""
"""
function plotVn()

end

