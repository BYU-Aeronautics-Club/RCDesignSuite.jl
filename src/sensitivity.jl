using PyPlot
using Roots
using FLOWMath


##########################################################################
####################         MISSION OBJECTIVES       ####################
##########################################################################

#=
This is where you input your mission objective functions.  Typically, the ground mission is hard to quantify without a lot of testing, but can be added.  In addition, the first flight mission is typically completion based, and therefore unnecessary to add here.
=#

"""
    mission2(inputs, normalization_factor=1.0)

Calculates score from mission 2 based on relavent inputs.

**Inputs:**
"""
function mission2(inputs=0.0, normalization_factor=1.0)

    # TODO: This is where you put the mission 2 function.  Note that the inputs variable is probably an array that you'll want to unpack before using.

    raw_score = inputs

    return raw_score/normalization_factor

end


"""
    mission3(inputs, normalization_factor=1.0)

Calculates score from mission 3 based on relavent inputs.

**Inputs:**
"""
function mission3(inputs=0.0, normalization_factor=1.0)

    # TODO: This is where you put the mission 3 function.  Note that the inputs variable is probably an array that you'll want to unpack before using.

    raw_score = inputs

    return raw_score/normalization_factor

end


"""
    ground_mission(inputs, normalization_factor=1.0)

!!OPTIONAL!!

Calculates score from ground mission based on relavent inputs.

**Inputs:**
"""
function ground_mission(inputs=1.0, normalization_factor=1.0)

    # TODO: This is where you put the mission 3 function.  Note that the inputs variable is probably an array that you'll want to unpack before using.

    raw_score = inputs

    return raw_score/normalization_factor

end



############################################################################
#####################        OBJECTIVE FUNCTION        #####################
############################################################################
"""
    objective(design_variables)

Function that takes in design variables and outputs the overall normalized mission score. (For the objective functions above.)
"""
function objective(design_variables; return_all=false)


    ### UNPACK DESIGN VARIABLES ###
    # TODO: Here is where you unpack the relevant design variables. For example:

    W               = design_variables[1] # Aircraft Weight
    S               = design_variables[2] # Wing Area
    CD0             = design_variables[3] # Zero Lift Drag Coeff of Aircraft
    eta             = design_variables[4] # Total Propulsive Efficiency
    P               = design_variables[5] # Available Battery Power

    #TODO: If you want to keep everything else the same, you can add more variables here.

    # Normalization Factors are the last 3 inputs in this example
    GM_norm_factor  = design_variables[end-2] # Ground Mission Normalization Factor
    M2_norm_factor  = design_variables[end-1] # Flight Mission 2 Normalization Factor
    M3_norm_factor  = design_variables[end] # Flight Mission 3 Normalization Factor

    ### GROUND MISSION ###
    # TODO: Here is where you set up the inputs for the ground mission if you desire.  Note that the function is set up to return 1.0 by default. (assumes perfect score)



    GM = ground_mission()


    ### FLIGHT MISSION 2 ###
    # TODO: Here is where you set up the inputs for flight mission 2 and call the objective function.



    M2 = mission2(M2_inputs, M2_norm_factor)


    ### FLIGHT MISSION 3 ###
    # TODO: Here is where you set up the inputs for flight mission 3 and call the objective function.



    M3 = mission3(M3_inputs, M2_norm_factor)



    ### SUM UP OBJECTIVES ###

    totalpoints = 100.0 + 1.0 + GM + M2 + M3 #assumes perfect score (100.0) on report and M1 completion (1.0)

    if return_all
        return totalpoints, GM, M2, M3
    else
        return totalpoints
    end

end



############################################################################
#####################        objective STUDY         #####################
############################################################################

# TODO: Define all the variables that might be important.

W               = 0.0 # Aircraft Weight
S               = 0.0 # Wing Area
CD0             = 0.0 # Zero Lift Drag Coeff of Aircraft
eta             = 0.0 # Total Propulsive Efficiency
P               = 0.0 # Available Battery Power



# TODO: Define input array for objective function.
design_variables        = zeros(8)
design_variables[1]     = W   # Aircraft Weight
design_variables[2]     = S   # Wing Area
design_variables[3]     = CD0 # Zero Lift Drag Coeff of Aircraft
design_variables[4]     = eta # Total Propulsive Efficiency
design_variables[5]     = P   # Available Battery Power
# Normalization Factors are the last 3 inputs in this example
design_variables[end-2] = 1.0 # Initialize GM_norm_factor to 1.0
design_variables[end-1] = 1.0 # Initialize M2_norm_factor to 1.0
design_variables[end]   = 1.0 # Initialize M3_norm_factor to 1.0



# Define range of percentages to do objective study over.
r = range(-50,stop=50,length=50)*1e-2

# Number of steps in the range
N = 50 #number of points in ranges


# Initialize the objective function outputs array
obj = zeros(length(design_variables),N)
# Initialize the objective function derivatives array
dobj = zeros(length(design_variables),N)


### GET NORMALIZATION FACTORS ###
# For each input (other than normalization factors), loop through defined range (r) from the nominal

# Initialize normalization factors.
GM_norm_factor = 0.0
M2_norm_factor = 0.0
M3_norm_factor = 0.0

# Run once to get maxima (or mean, or whatever you want to normalize by) for mission normalization factors.
for i = 1:length(design_variables)-3 #(Don't include norm factors)

    # Create a copy of the original design variables so you can isolate each one in turn.
    design_variablescopy = copy(design_variables)

    for j = 1:N

        # Loop through each of the design variables across the defined relative range.
        design_variablescopy[i] = (1.0 + r[j])*design_variables[i]

        #try catch in case you get a combination that doesn't make sense or has a negative square root or something.
        try

            #throw away the total, we don't need it, we're just getting the maximum components for normalization now.
            _, GM_temp, M2_temp, M3_temp = objective(design_variablescopy)

            # Find Maximum normalization values.
            if GM_temp > GM_norm_factor
                GM_norm_factor = GM_temp
            end
            if M2_temp > M2_norm_factor
                M2_norm_factor = M2_temp
            end
            if M3_temp > M3_norm_factor
                M3_norm_factor = M3_temp
            end

        catch

            #printl the cases that don't work so you can debug and make adjustments as needed.
            println("i: $i, j: $j")

        end #try/catch
    end #for range
end #for number of design variables


## Put your Normalization Factors in your Design Variables Array:
# Normalization Factors are the last 3 inputs in this example
design_variables[end-2] = GM_norm_factor
design_variables[end-1] = M2_norm_factor
design_variables[end] = M3_norm_factor



### GET SENSITIVITIES ###

# Run a second time to get actual values for the objectives, normalized by the max GM, M2, and M3 values just found.
for i = 1:length(design_variables)-3 #(Don't include norm factors)

    # Same copying as before
    design_variablescopy = copy(design_variables)

    for j = 1:N

        # Same loop over range as before.
        design_variablescopy[i] = (1.0 + r[j])*design_variables[i]

        try

            # Don't return everything this time.
            obj[i,j] = objective(design_variablescopy)

        catch

            # Set failures to NaN so they don't show up in the plot.
            obj[i,j] = NaN

        end #try/catch
    end #for range

    # Spline the outputs of each design variable range.
    so = FLOWMath.Akima(r, obj[i,:])

    # Get the derivatives relative to the percent of input.
    for j = 1:N
        dobj[i,j] = FLOWMath.derivative(so,r[j])
    end

end # for number of variables


## Calculate the objective value for the nominal case for getting the percent differences later.
obj0 = objective(design_variables)



############################################################################
#####################       INTERMEDIATE PLOTS         #####################
############################################################################

# TODO: Create a vector for the labels that match your design variables, so you know what the following plots are for.
labels = [
"W";
"S";
"CD0";
"eta";
"P";]

# Plot each variable on it's own plot so you can clearly see its sensitivity.
for i=1:length(design_variables)-3
    figure()
    clf()
    plot(r,(obj[i,:].-obj0)./obj0,label=labels[i])
    legend()
end



############################################################################
#####################           FINAL PLOTS            #####################
############################################################################

#=
HERE YOU NEED TO FIRST LOOK AT ALL THE INDIVIDUAL PLOTS AND SEE WHICH ONES AFFECT THE OBJECTIVE
THEN YOU NEED TO SEE IF ANY ARE THE EXACT SAME (NO NEED TO PLOT LINES ATOP EACHOTHER)
THEN YOU NEED TO PICK WHICH ONES YOU WANT ON THE FINAL PLOT
=#

# TODO: identify the index of the variables you want
toplot = [1;2;3;4;5]

# Get the plot labels you want to show for those variables
toplotlabels = labels[toplot]

# Pick your colors
colors = ["C0";"C1";"C2";"C0";"C1";"C2"]

# Pick your styles.
styles = ["-";"-";"-";"--";"--";"--"]

# If needed, add a scale factor that helps the plot be more clearly
scalefactor = 1e3



#### PLOT FINAL FIGURE FOR REPORTS ####
# Initialize the figure so that it will be about the right size for the paper (4x3 ratio), but also such that the legend doesn't overlap things weird.
figure(figsize=(4.5,4.5*3/4))

xlabel("Percent Change in Design Variable (from nominal)")
ylabel("Normalized Change in Score (x $scalefactor)")

# Plot the interesting things.
for i=1:length(toplot)
    plot(r*1e2,scalefactor*(obj[toplot[i],:].-obj0)./obj0,label=toplotlabels[i],linestyle=styles[i],color=colors[i])
end

legend()

#save the figure.
savefig("../figs/sensitivityobj.png",bbox_inches="tight")



#### MAYBE PLOT THE DERIVATIVES INSTEAD, IF THAT MAKES MORE SENSE ####
#initialize the figure so that it will be about the right size for the paper (4x3 ratio), but also such that the legend doesn't overlap things weird.
figure(figsize=(4.5,4.5*3/4))

xlabel("Percent Change in Design Variable (from nominal)")
ylabel("Differential Change in Score (dScore/dx)")

#plot the interesting things.
for i=1:length(toplot)
    plot(r*1e2,dobj[toplot[i],:],label=toplotlabels[i],linestyle=styles[i],color=colors[i])
end

legend()

#save the figure.
savefig("../figs/sensitivitydobj.png",bbox_inches="tight")
