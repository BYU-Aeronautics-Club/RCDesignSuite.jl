#=

Authors: Judd Mehr,

Preliminary Design functions

=#

export addairfoilflaps, saveairfoilpolar, airfoilcomp


############################
#######    General    ######
############################



############################
#######    Sizing    #######
############################



############################
####    Aerodynamics    ####
############################
"""
    addairfoilflaps(filename,datapath,xc,yt,angle;npan=140)

Add a flap to airfoil geometry at a given percent chord and percent thickness at that location.

*Inputs:*
- filename::String : file name of airfoil geometry (format: column 1 = x/c, column 2 = y/c)
- suffix::String : file type of file (e.g. .dat, .csv, etc.)
- datapath::String : the path where the file can be found
- xc::Float64 : x/c location where flap is to be defined (e.g. 0.7 is at 70% chord)
- yt::Float64 : y/t relative to thickness at x/c location (e.g. 0.5 is halfway between upper and lower surfaces).
- angle::Float64 : angle of flap deflection, in degrees. + is down, - is up.

*Keyword Arguments:*
- npan::Int64 : number of panels to use for geometry
- savefile::String : path and file name where to save output geometry.

*Outputs:*
- geometry file with flap added, default file name is datapath*filename*"$(angle)deg_$(xc*100)chord.dat"
"""
function addairfoilflaps(filename, suffix, datapath, xc, yt, angle;
                npan=140,
                savefile=nothing
                )

    #Get original airfoil data
    xy = convert(Array,DataFrames.DataFrame!(CSV.File(datapath*filename*suffix,header=false,datarow=2)))

    #Pane data so it's standardized
    Xfoil.set_coordinates(xy[:,1],xy[:,2])
    Xfoil.pane(npan=npan)
    x = Xfoil.get_globals().x[1:npan]
    y = Xfoil.get_globals().y[1:npan]

    #identify points to rotate
    flapidx = findall(x->x>=xc,x)
    flaparr = [x[flapidx] y[flapidx]]

    #calculate y-axis rotation point
    #find halfway point of flapidx
    halfway = Int(length(flapidx)/2)

    #then grab y points there and find distance between
    ylo = y[flapidx[halfway]]
    yhi = y[flapidx[halfway+1]]
    #then calculate actual y point based on relative thickness input
    ypt = ylo + yt*(yhi-ylo)

    #translate to origin from rotation point
    flaparr .-= [xc ypt]

    #rotate
    rotmat = [cosd(-angle) -sind(-angle);
              sind(-angle)  cosd(-angle)]

    rotarr = rotmat*flaparr'
    rotarr = rotarr'

    #translate back
    rotarr .+= [xc ypt]

    #replace previous flap locations with rotated ones
    xnew = copy(x)
    xnew[flapidx] = rotarr[:,1]
    ynew = copy(y)
    ynew[flapidx] = rotarr[:,2]

    #delete the points that dip below (or above) old y values to avoid kinks.
    if sign(angle) < 0 #if flap is up
        rmidx = findall(x->x<y[flapidx[halfway]],ynew[flapidx[1:halfway]])
        deleteat!(ynew,flapidx[rmidx])
        deleteat!(xnew,flapidx[rmidx])
    else
        rmidx = findall(x->x>y[flapidx[halfway+1]],ynew[flapidx[halfway+1]:end])
        deleteat!(ynew,flapidx[halfway.+rmidx])
        deleteat!(xnew,flapidx[halfway.+rmidx])
    end

    #savefile:
    if savefile === nothing
        savefile =  datapath*filename*"$(angle)deg_$(xc*100)chord.dat"
    end

    f = open(savefile, "w")

    string = "x/c y/c\n"

    write(f, string)

    for i=1:length(xnew)

        string = "$(xnew[i]) $(ynew[i])\n"

        write(f, string)

    end

    close(f)

end


"""
    saveairfoilpolar(filename, datapath, Re, alpha, cl, cd, cdp, cm, conv)

Save outputs from Xfoil.jl analysis.

*Inputs:*
- filename::String : file name of save file.
- datapath::String : path where to save file.
- Re::Float64 : Reynolds number.
- alpha::Array{Float64} : angles of attack.
- cl::Array{Float64} : lift coefficients.
- cd::Array{Float64} : drag coefficients.
- cdp::Array{Float64} : parasitic drag coefficients.
- cm::Array{Float64} : moment coefficients.
- conv::Array{Bool} : convergence flags.

*Keyword Arguments:*
- savefile::String : path and file where to save output.

*Outputs:*
- saved file with location and name. Default: datapath*filename*"$(@sprintf("%.2E", Re)).csv"  File will always be saved with "$(@sprintf("%.2E", Re)).csv" tagged on the end.
"""
function saveairfoilpolar(filename, datapath, Re, alpha, cl, cd, cdp, cm, conv;
                    savefile = nothing)

    if savefile === nothing
        savefile = datapath*filename*"$(@sprintf("%.2E", Re)).csv"
    else
        savefile *= "$(@sprintf("%.2E", Re)).csv"
    end

    f = open(savefile, "w")

    string = "alpha,cl,cd,cdp,cm,conv\n"

    write(f, string)

    for i=1:length(alpha)

        string = "$(alpha[i]),$(cl[i]),$(cd[i]),$(cdp[i]),$(cm[i]),$(conv[i])\n"

        write(f, string)

    end

    close(f)

end

"""
    runxfoil(filenames, datapath, Revalues, alpharange=collect(-5:20); savefile)

Run xfoil for give airfoil geometry and Reynolds numbers, and pass data to saveairfoilpolar().

*Inputs:*
- filenames::Array{String} : geometry file names of airfoils to analyze.
- datapath::String : path where files are located.
- alpharange::Array{Float64} : angles of attack (in degrees) to analyze.

*Keyword Arguments:*
- savefile::Array{String} : names of files (including path) to save the polar data. (File will always be saved with "$(@sprintf("%.2E", Re)).csv" tagged on the end.)
"""
function runxfoil(filenames, datapath, Revalues, alpharange=collect(-5:20);
            savefile = nothing)

    for i=1:length(filenames)

        #Load Airfoil
        xy = convert(Array,DataFrames.DataFrame!(CSV.File(datapath*filenames[i],header=false,datarow=2)))

        #Loop through Reynolds numbers
        for j=1:length(Revalues)

            #Run Xfoil
            cl, cd, cdp, cm, conv = Xfoil.alpha_sweep(xy[:,1],xy[:,2], alpharange, Revalues[j])

            #Save data so you don't have to re-run xfoil
            saveairfoilpolar(filenames[i],datapath,Revalues[j],alpharange,cl,cd,cdp,cm,conv)

        end

    end

end


############################
#####    Propulsion    #####
############################



############################
#####    Structures    #####
############################