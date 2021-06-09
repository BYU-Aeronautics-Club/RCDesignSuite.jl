#=
Authors: Judd Mehr,

Several ready-made plotting functions for common plots
=#


export standardizeafs, plotaf, plotpolar, ploteigenvals, plotbasic, plotpropefficiency, plotdragbars, plotliftdist, plotLoverDvVel, plotLoverDvCL, plotLoverD



#######################
####### General #######
#######################

"""
Plots the induced, parasitic, and total drag with repsect to velocity.
Inputs
* A struct of type `concept`
"""
function dragplot(a::concept)
    S = a.area
    w = a.weight
    v, Dp, Di, Dt = dragdata(a)
    plot(v,Dp, label = "Parasitic", xlabel = "Velocity (m/s)", ylabel = "Drag (N)")
    plot!(v,Di, label = "Induced")
    plot!(v,Dt, label = "Total")
end

"""
Plots lift coefficient with repsect to velocity.
Inputs
* A struct of type `concept`
"""
function plotliftcoeff(a::concept)
    S = a.area
    w = a.weight
    v, CL, _ = liftcoefficient(a)
    plot(v,CL, label = "CL", xlabel = "Velocity (m/s)", ylabel = "CL")
end

"""
    readxflrgraphs(filename)

Reads in arrays from XFLR5 graph export .csv files.

*Inputs:*
- filename::String : location and file name of data file.

*Outputs:*
- xs::Array{Float64} : "x" values of input data
- ys::Array{Float64} : "y" values of input data
"""
function readxflrgraphs(filename)

    #get the real and imaginary parts
    mat = convert(Array,DataFrames.DataFrame!(CSV.File(filename,header=false,datarow=2,delim=",")))

    #get raw x components
    xs = mat[:,1:3:end]

    #remove last column of missing's
    if typeof(xs[end])==Missing
        xs = xs[:,1:end-1]
    end

    #get raw y components
    ys   = mat[:,2:3:end]

    #clear out other missings
    xs = replace(xs,missing=>NaN)
    ys = replace(ys,missing=>NaN)

    return xs, ys

end



"""
    plotbasic(xs, ys, labels)

Create basic plot with default settings

*Inputs:*
- xs::Array{Float64} : horizontal axis values
- ys::Array{Float64} : vertical axis values
- labels::Array{String} : labels for each dataset

*Keyword Arguments:*
- savefigure::Bool : save figure if true
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (4,3)
"""
function plotbasic(xs, ys, labels;
            savefigure  = true,
            savefile    = nothing,
            savepath    = nothing,
            wh          = (4,3)
            )

    plt.figure(figsize=wh)
    for i=1:length(xs[1,:])
        plt.plot(xs[:,i],ys[:,i],labels[i])
    end

    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "basic_plot.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

    return false
end



"""
for plotting against testing data...
"""
function ploterrbars()

end



#######################
###### AIRFOILS #######
#######################
"""
    standardizeafs(filenames, datapath)

Imports airfoil files and makes the data arrays the same length for correct inputs to plotaf.

*Inputs:*
- filenames::Array{String} : Array of the airfoil filenames
- datapath::String : String for path where airfoil files can be found

*Outputs:*
- xs::Array{Float64} : Array of x coordinates of the airfoil data, each column is an airfoil.
- ys::Array{Float64} : Array of y coordinates of the airfoil data, each column is an airfoil.
"""
function standardizeafs(filenames, datapath)

    #Read in first airfoil corrdinates
    mat = DataFrames.DataFrame!(CSV.File(datapath*filenames[1],header=true,datarow=2))

    #get array lengths
    ldata = length(mat[:,1])
    lafs = length(filenames)

    #Initialize x and y arrays
    xs = zeros(ldata,lafs)
    ys = zeros(ldata,lafs)

    #fill in first columns
    xs[:,1] = mat[:,1]
    ys[:,1] = mat[:,2]

    #Loop through the rest of the files
    for i=2:lafs

        #get current length of xs
        ldata = length(xs[:,1])

        #Read in next airfoil coordinates
        tempmat = DataFrames.DataFrame!(CSV.File(datapath*filenames[i],header=true,datarow=2))
        ltemp = length(tempmat[:,1])
        tempx = tempmat[:,1]
        tempy = tempmat[:,2]

        #If the current airfoil data is shorter than the rest, fill the gap with NaNs (that won't plot)
        if ltemp < ldata
            tempx = [tempx; NaN*ones(ldata-ltemp)]
            tempy = [tempy; NaN*ones(ldata-ltemp)]

        #If the current airfoil data is longer, fill the gap with NaNs (that won't plot)
        elseif ltemp > ldata
            xs = [xs; NaN*ones(ltemp-ldata,lafs)]
            ys = [ys; NaN.*ones(ltemp-ldata,lafs)]
        end

        #Add the current data into the full set.
        xs[:,i] = tempx
        ys[:,i] = tempy

    end

    return xs, ys

end



"""
    plotaf(x, y, labels)

Plot airfoil(s) with unit chord.

Note: the airfoil coordinate vectors need to be the same to concatenate them. Simply add the final point repeatedly to the shorter vector until  long enough.

*Inputs:*
- x::Array{Float64} : Array of x coordinates for the airfoil(s). (each airfoil is a column of x)
- y::Array{Float64} : Array of y coordinates for the airfoil(s). (each airfoil is a column of y)
- labels::Array{String} : Array of labels (airfoil names)

*Keyword Arguments:*
- savefigure::Bool : boolean of whether or not to save the figure automatically.
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- width::Float64 : width of figure (in inches), default = 7
- colorcycle::Array{String} : array of colors for plot. Must be >= number of airfoils.
- stylecycle::Array{String} : array of line styles for plot. Must be >= number of airfoils.

*Outputs:*
- plot of whatever the input width is and height proportional to 5% greater than the relative height of the maximum y-location of all the airfoils.
"""
function plotaf(filenames, datapath, labels;
            savefigure  = true,
            savefile    = nothing,
            savepath    = nothing,
            width       = 7,
            colorcycle  = ["C0"; "C1"; "C3"; "C0"; "C1"; "C3"], stylecycle  = ["-";"-";"-";"--";"--";"--"]
            )

    x, y = standardizeafs(filenames, datapath)

    if length(x[1,:]) > length(colorcycle) || length(x[1,:]) > length(stylecycle)
        @error("You are trying to plot more than 6 airfoils, you need to define a longer color and/or style cycle.")
    end

    maxx = maximum(filter(!isnan,x))
    maxy = maximum(filter(!isnan,y))

    height = width*2.1*maxy
    plt.figure(figsize=(width,height))
    for i=1:length(x[1,:])
        maxxs = maximum(filter(!isnan,x[:,i]))
        plt.plot(x[:,i]./maxxs,y[:,i]./maxxs,color=colorcycle[i],linestyle=stylecycle[i],label=labels[i])
    end
    plt.legend()
    plt.axis(false)

    #save figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "airfoil_geometry.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

    return false

end



"""
    getticks(data, spacing)

Find plotting ticks based on data range and chosen spacing.

*Inputs:*
- data::Array{Float64} : data array to define plotting ticks for.
- spacing::Float64 : value of tick spacing.

*Outputs:*
- ticks::Array{Float64} : Array of plot ticks.
"""
function getticks(data, spacing)

    tickloc = []
    for i=1:length(data)
        if data[i] != NaN
            if data[i] == 0.0
                push!(tickloc,0.0)
            else
                push!(tickloc,sign(data[i])*round(Int, abs(data[i]) / spacing) * spacing)
            end
        end
    end

    return unique(tickloc)

end



"""
    plotafpolar(clvalphafile, clvcdfile, datapath, labels)

Reads in XFLR5 airfoil polar plot files and plots them nicely.

Assumes you're using polars from all the same XFLR5 project (otherwise things won't match up).

*Inputs:*
- clvalphafile::String : Airfoil polar file name for cl vs alpha
- cdvclfile::String : Airfoil polar file name for cl vs cd
- datapath::String : String for path where airfoil polar files can be found
- labels::Array{String} : Array of labels (airfoil names)

*Keyword Arguments:*
- savefigure::Bool : boolean of whether or not to save the figure automatically.
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (4,3)
"""
function plotpolar(clvalphafile, clvcdfile, datapath, labels;
                savefigure  = true,
                savefile    = nothing,
                savepath    = nothing,
                wh          = (7,2.75),
                alphaticspace= 3.0,
                cdticspace  = 0.03,
                airfoil     = true,
                )

    #Initialize Figure
    plt.figure(figsize=wh)

    #c_l vs alpha subplot
    alpha, cl = readxflrgraphs(datapath*clvalphafile)
    ax1 = plt.subplot(131)
    plt.xlabel(plt.L"Angle of Attack ($\alpha^\circ$)")
    if airfoil
        plt.ylabel(plt.L"Lift Coefficient ($c_\ell$)")
    else
        plt.ylabel(plt.L"Lift Coefficient ($C_L$)")
    end

    for i=1:length(alpha[1,:])
        plt.plot(alpha[:,i],cl[:,i],label=labels[i])
    end
    plt.plot([0.0;0.0],[minimum(cl),maximum(cl)],":C2")
    xtickloc = getticks(alpha, alphaticspace)
    ax1.set_xticks(xtickloc)


    #c_l vs c_d subplot
    cd, cl = readxflrgraphs(datapath*clvcdfile)
    ax2 = plt.subplot(132)

    if airfoil
        plt.ylabel(plt.L"Lift Coefficient ($c_\ell$)")
        plt.xlabel(plt.L"Drag Coefficient ($c_d$)")
    else
        plt.ylabel(plt.L"Lift Coefficient ($C_L$)")
        plt.xlabel(plt.L"Drag Coefficient ($C_D$)")
    end

    for i=1:length(cd[1,:])
        plt.plot(cd[:,i],cl[:,i],label=labels[i])
    end
    # plt.plot([0.0;0.0],[minimum(cl),maximum(cl)],":C2")
    xtickloc = getticks(cd, cdticspace)
    ax2.set_xticks(xtickloc)


    #legend subplot
    plt.subplot(133)
    for i=1:length(alpha[1,:])
        plt.plot([NaN],[NaN],label=labels[i])
    end
    plt.axis(false)
    plt.legend(loc="center left")

    #save figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "polar.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

end



#######################
###### AIRCRAFT #######
#######################
"""
    plotliftdist(filename, datapath, totalspan, labels; args...)

Plot wing lift distribution.

*Inputs:*
*Inputs:*
- filename::String : Filename of XFLR5 export file.
- datapath::String : String for path where airfoil file can be found
- totalspan::Float64 : The total span value of the wing, in the same units as the datafile.
- labels::Array{String} : array of labels for each of the lift distributions contained in the datafile.


*Keyword Arguments:*
- savefigure::Bool : boolean of whether or not to save the figure automatically.
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (6,3)
"""
function plotliftdist(filename, datapath, totalspan, labels;
            savefigure  = true,
            savefile    = nothing,
            savepath    = nothing,
            wh          = (6,3))

    #Read data
    spanloc, cl = readxflrgraphs(datapath*filename)

    spannorm = (spanloc .+ abs(minimum(spanloc))) ./ totalspan

    plt.figure(figsize=wh)
    plt.ylabel(plt.L"Local Lift Coefficient ($c_\ell$)")
    plt.xlabel(plt.L"Normalized Spanwise Location ($x/b$)")

    for i=1:length(spanloc[1,:])
        plt.plot(spannorm,cl,label=labels[i])
    end
    plt.legend()

    #save the figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "lift_distribution.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end
end



"""
    plotLoverDvVel(filename,datapath,labels; args...)

Plot L/D vs Airspeed.

*Inputs:*
- filename::String : name of XFLR export file containing L/D vs V data (L/D on y-axis, V on x-axis)
- datapath::String : path where file is saved.
- labels::Array{String} : Array of configuration names in data file.

*Keyword Arguments:*
- designvelocity::Float64 : Design Velocity value (plots as vertical  line).
- savefigure::Bool : boolean of whether or not to save the figure automatically.
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (4,3)
"""
function plotLoverDvVel(filename,labels;
    datapath        = "./data/aircraft/",
    designvelocity  = nothing,
    units           = "metric",
    savefigure      = true,
    savefile        = nothing,
    savepath        = nothing,
    wh              = (4,3))

    #Read data
    V, LD = RCDesignSuite.readxflrgraphs(datapath*filename)

    if units == "metric"
        V *= 3.28084 #convert to ft/s from m/s
        if designvelocity != nothing
            designvelocity *= 3.28084
        end
    end

    plt.figure(figsize=wh)
    plt.ylabel("Lift to Drag Ratio (L/D)")
    plt.xlabel("Velocity (ft/s)")

    for i=1:length(V[1,:])
        plt.plot(V[:,i],LD[:,i],label=labels[i])
    end
    if designvelocity !== nothing
        plt.plot([designvelocity; designvelocity], [minimum(LD); maximum(LD)],"--",label="Cruise Velocity")
    end
    plt.legend()

    #save the figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "L_over_D_vs_Vel.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

end



"""
    plotLoverDvCL(filename,datapath,labels; args...)

Plot L/D vs CL.

*Inputs:*
- filename::String : name of XFLR export file containing CL vs CD data (CD on x-axis, CL on y-axis)
- datapath::String : path where file is saved.
- labels::Array{String} : Array of configuration names in data file.

*Keyword Arguments:*
- designcl::Float64 : Design CL value (plots as vertical  line).
- savefigure::Bool : boolean of whether or not to save the figure automatically.
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (4,3)
"""
function plotLoverDvCL(filename,datapath,labels;
    designcl = nothing,
    savefigure  = true,
    savefile    = nothing,
    savepath    = nothing,
    wh          = (4,3))

    #Read data
    CD, CL = readxflrgraphs(datapath*filename)
    LD = CL./CD

    plt.figure(figsize=wh)
    plt.ylabel("Lift to Drag Ratio (L/D)")
    plt.xlabel(plt.L"Lift Coefficient ($C_L$)")

    for i=1:length(CL[1,:])
        plt.plot(CL[:,i],LD[:,i],label=labels[i])
    end
    if designcl !== nothing
        plt.plot([designcl; designcl], [minimum(LD); maximum(LD)],"--",label="Design Lift Coefficient")
    end
    plt.legend()

    #save the figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "L_over_D_vs_CL.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

end

"""
"""
function plotLoverD(filename,labels;
    datapath    = "./data/aircraft/",
    savefigure  = true,
    savefile    = nothing,
    savepath    = nothing,
    wh          = (4,3))

    #Read data
    CD, CL = readxflrgraphs(datapath*filename)

    plt.figure(figsize=wh)
    plt.xlabel(plt.L"Drag Coefficient ($C_D$)")
    plt.ylabel(plt.L"Lift Coefficient ($C_L$)")

    for i=1:length(CL[1,:])
        plt.plot(CD[:,i],CL[:,i],label=labels[i])
    end
    plt.legend()

    #save the figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "CL_v_CD.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

end

#######################
###### STABILITY ######
#######################
"""
    removenegatives!(checkarray,followarray)

Checks for negative values in the first array and deletes the idices from both arrays.

*Inputs:*
- checkarray::Array : the array to check for negative values
- followarray::Array : the other array to remove entries from

*Outputs:*
- checkarray::Array : the array to check for negative values, with negatives removed
- followarray::Array : the other array to remove entries from, with matching indices removed
"""
function removenegatives!(checkarray,followarray)
    #get indices for negative ys
    negidx = findall(x->x<0,checkarray)

    #remove negative imaginary entries
    deleteat!(checkarray,negidx)
    deleteat!(followarray,negidx)

end



"""
    ploteigenvals(latfilename, lonfilename, datapath; args...)

Read csv files from XFLR5 exported eigenvalue plots and plot the eigenvalues nicely.

It is assumed that all configurations are contained in the same file (as what happens when you plot the various configurations on a single plot in XFLR5)

*Inputs:*
- lonfilename::Array{String} : array of longitudinal file names
- latfilename::Array{String} : array of lateral file names
- datapath::String : path to where eigenvalue files are saved

*Keyword Inputs:*
- savefigure::Bool : save figure if true
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (6,4)
- lonmodenames::Array{String} : array of longitudinal mode names for plot annotation, default = ["Short Period";"Phugoid"]
- latmodenames::Array{String} : array of lateral mode names for plot annotation, default = ["Roll";"Dutch Roll";"Spiral"]
"""
function ploteigenvals(nconfig, latfilename, lonfilename, labels;
                datapath        = "./data/stability/",
                savefigure      = true,
                savefile        = nothing,
                savepath        = nothing,
                wh              = (6,4),
                lonmodenames    = ["Short Period";"Phugoid"],
                latmodenames    = ["Roll";"Dutch Roll";"Spiral"],
                colorcycle      = ["C0"; "C2"; "C3"],
                markercycle     = ["s"; "^"; "."]
                )

    latr = zeros(3,nconfig)
    lati = zeros(3,nconfig)
    lonr = zeros(2,nconfig)
    loni = zeros(2,nconfig)

    #Read and define the 3 lateral modes
    res, ims = RCDesignSuite.readxflrgraphs(datapath*latfilename)
    #just grab the first row, since everything else is repeats
    res = res[1,:]
    ims = ims[1,:]
    #remove negatives so we just plot one quadrant (since it's symmetric)
    RCDesignSuite.removenegatives!(ims,res)

    #load arrays for plotting
    for i=1:nconfig
        latr[:,i], lati[:,i] = res[3*(i-1)+1:3*i], ims[3*(i-1)+1:3*i]
    end

    #Repeat for the 2 logitudinal modes
    res, ims = RCDesignSuite.readxflrgraphs(datapath*lonfilename)
    res = res[1,:]
    ims = ims[1,:]
    RCDesignSuite.removenegatives!(ims,res)
    for i=1:nconfig
        lonr[:,i], loni[:,i] = res[2*(i-1)+1:2*i], ims[2*(i-1)+1:2*i]
    end

    #initialize plot
    plt.figure(figsize=wh)
    plt.xlabel(plt.L"Re($\lambda$)")
    plt.ylabel(plt.L"Im($\lambda$)")

    #lateral modes
    for j=1:3
        for i=1:nconfig
            # if latr[j,i] > 0
            #     color = "C1"
            # else
            color = colorcycle[i]
            marker = markercycle[i]
            # end

            if j==1
                plt.plot(latr[j,i],lati[j,i],color=color,marker=marker,label=labels[i])
            else
                plt.plot(latr[j,i],lati[j,i],color=color,marker=marker,label="")
            end
        end
        plt.annotate(latmodenames[j],xy=[maximum(latr[j,:])+0.5,maximum(lati[j,:])-0.03125],color="k")
    end

    #lonitudinal modes
    for j=1:2
        for i=1:nconfig
            # if lonr[j,i] > 0
            #     color = "C1"
            # else
            color = colorcycle[i]
            marker = markercycle[i]
            # end

            if j==1
                plt.plot(lonr[j,i],loni[j,i],color=color,marker=marker,label="")
            else
                plt.plot(lonr[j,i],loni[j,i],color=color,marker=marker,label="")
            end
        end
        plt.annotate(lonmodenames[j],xy=[maximum(lonr[j,:])+0.5,maximum(loni[j,:])-0.03125],color="k")
    end

    #set axis bounds
    ylow = -0.25
    yhi = max(maximum(lati),maximum(loni))+0.5
    plt.ylim([ylow,yhi])
    xlow = min(minimum(latr),minimum(lonr))-0.5
    xhi = 4.9
    plt.xlim([xlow,xhi])

    #add a dotted line at the boundary of the left half plane
    plt.plot([0.0;0.0],[ylow;yhi],":C2",zorder=0)

    #move the axes over to the right so it makes more sense to look at.
    ax = plt.gca()
    ax[:spines]["left"][:set_visible](false)
    ax[:spines]["right"][:set_visible](true)
    ax[:yaxis][:set_label_position]("right")
    ax[:yaxis][:set_ticks_position]("right")
    plt.legend()

    #save the figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "eigen_values.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

end



"""
    plotcm(filename, datapath, labels; args...)

Reads in XFLR5 cm vs alpha plot file and plots it nicely.

*Inputs:*
- filename::String : Aircraft polar file name for Cm vs alpha
- datapath::String : String for path where polar file can be found
- labels::Array{String} : Array of labels of configurations contained in data file

*Keyword Arguments:*
- savefigure::Bool : boolean of whether or not to save the figure automatically.
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (4,3)
"""
function plotcm(filename,datapath,labels;
    savefigure  = true,
    savefile    = nothing,
    savepath    = nothing,
    wh          = (4,3))

    #Read data
    alpha, cm = readxflrgraphs(datapath*filename)

    plt.figure(figsize=wh)
    plt.ylabel(plt.L"Moment Coefficient ($C_M$)")
    plt.xlabel(plt.L"Angle of Attack ($\alpha^\circ$)")

    for i=1:length(alpha[1,:])
        plt.plot(alpha,cm,label=labels[i])
    end
    plt.legend()

    #save the figure
    if savefigure

        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "lift_distribution.pdf"
        end

    plt.savefig(savepath*savefile, bbox_inches="tight")

    end
end



#######################
##### PERFORMANCE #####
#######################
"""
    plotpropefficiency(J,etamotor,etaprop)

Plot propulsion efficiency.

*Inputs:*
- J::Array{Float64} : array of advance ratios
- etamotor::Array{Float64} : array of motor efficiencies
- etaprop::Array{Float64} : array of propeller efficiencies

*Keyword Inputs:*
- savefigure::Bool : save figure if true
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (6,4)
"""
function plotpropefficiency(J,etamotor,etaprop;
                    savefigure = true,
                    savefile   = nothing,
                    savepath   = nothing,
                    wh         = (4,3)
                    )

    plt.figure(figsize=wh)
    plt.plot(J,etamotor,label="Motor",color="C0")
    plt.plot(J,etaprop,label="Prop",color="C2")
    plt.plot(J,etamotor.*eteaprop,label="Total",color="C3")
    plt.xlabel("Advance Ratio (J)")
    plt.ylabel(plt.L"Efficiency ($\eta$)")

    #save the figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "prop_efficiency.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

end



"""
    plotdragbars(dragcomp, labels)

Plot drag components as horizontal stack bars.

*Inputs:*
- dragcomp::Array{Float64,2} : array of drag componenets for each mission in percent.
- labels::Array{Float64} : array of labels for the drag components

*Keyword Inputs:*
- savefigure::Bool : save figure if true
- savefile::String : filename of figure to save
- savepath::String : path to location of where to save file
- wh::Tuple : width and height of figure, default = (5,2)
- width::Float64 : width of bars
"""
function plotdragbars(dragcomp, labels;
                savefigure  = true,
                savefile    = nothing,
                savepath    = nothing,
                wh          = (5,2),
                width       = 0.5
                )

    ind = 1:3   # the x locations for the groups

    plt.figure(figsize=wh)
    for i=1:length(dragcomp[1,:])
        plt.barh(ind, dragcomp[:,i], width,label=labels[i])
    end
    plt.xlabel("Percent Drag")
    plt.yticks(ind, ("M1", "M2", "M3"))
    plt.xticks(0:10:100)
    plt.legend()

    #save the figure
    if savefigure
        if savepath === nothing
            savepath = "./figs/"
        end

        if savefile === nothing
            savefile = "prop_efficiency.pdf"
        end

        plt.savefig(savepath*savefile, bbox_inches="tight")

    end

end



"""
"""
function plotVn()

end

