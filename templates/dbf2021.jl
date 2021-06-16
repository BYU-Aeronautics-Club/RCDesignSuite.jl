#=

Previous Years setup, objective, and constraint functions for examples and comparisons.

Authors: Judd Mehr,

=#

# using PyPlot, RCDesignSuite

#############################################
###########          2021         ###########
#############################################
function setup2021()

    # Initial Values
    x0 = [
        25.0/10.0    # Cruise Velocity, m/s
        10.0/10.0    # number of containers
        0.25    # sensor length, meters
        0.25    # individual sensor weight, kg
        1.0     # battery mass, kg
        25.0/10.0    # battery C rating
        21.0/10.0    # battery voltage
        0.25    # wing area, m^2
    ]

    # Lower Bounds
    lb = [
        0.0/10.0     # Cruise Velocity, m/s
        1.0/10.0     # number of containers
        0.005   # sensor length, meters
        0.0     # individual sensor weight, kg
        0.0     # battery mass, kg
        10.0/10.0    # battery C rating
        3.7/10.0     # battery voltage
        0.0     # wing area , m^2
    ]

    # Upper Bounds
    ub = [
        30.0/10.0   # Cruise Velocity, m/s
        20.0/10.0   # number of containers
        0.5     # sensor length, meters
        1.0     # individual sensor weight, kg
        1.0     # battery mass, kg
        50.0/10.0    # battery C rating
        34.0/10.0    # battery voltage
        2.0     # wing area, m^2
    ]

    # Parameter Values
    p = [
        1200.0  # assumed length of lap, meters
        1.0     # normalization factor for mission 2
        1.0     # normalization factor for mission 3
        75.0    # battery specific energy, watt-hours/kg
        0.6     # propulsive efficiency
        0.8     # wing CL max
        1.225   # air density, kg/m^3
        0.50    # M2 empty weight ratio (ratio of empty weight to total weight)
        0.75    # M3 empty weight ratio
        9.81    # gravity acceleration constant
        10.0    # lift to drag ratio
        0.05    # zero lift drag coefficient
    ]

    # Constraint Values
    c = [
        10.0   # max time, in minutes, allowed for mission 3 flight.
        1.524   # max wingspan, meters (projected)
        20.0    # maximum allowed weight, kg
        30.0    # maximum takeoff distance, meters
        200.0   # maximum battery stored power watt-hours
        0.10    # maximum  root bending moment (need to find an actaul number here.)
        20.0    # maximum aspect ratio
    ]


    _, m2, m3 = obj2021(ub, p, c; return_all=true)

    p[2] = m2
    p[3] = m3

    return x0, lb, ub, p, c

end



function obj2021(x, p, c; return_all=false)

    ### --- Unpack Variables
    # Unpack applicable design variables
    cruisevelocity  = x[1]*10.0 # cruise velocity, feet/sec
    ncontainers     = x[2]*10.0 # number of containers
    sensorlength    = x[3] # sensor length
    sensorweight    = x[4] # sensor weight

    # Unpack applicable parameters
    laplength       = p[1] # assumed length of lap
    M2_norm_factor  = p[2] # normalization factor for mission 2
    M3_norm_factor  = p[3] # normalization factor for mission 3

    # Unpack applicable constraints
    ttotal          = c[1] # max time allowed for mission 3



    ### --- Calculate intermediate values
    # time per lap
    tperlap = laplength/cruisevelocity

    # time for 3 laps on mission 2
    t2 = 3*tperlap

    # number of total  laps for mission 3
    nlaps = ttotal/tperlap



    ### --- Compile Objectives
    # Flight Mission 2 Score
    M2 = 1.0 + (ncontainers/t2) / M2_norm_factor

    # Flight Mission 3 Score
    M3 = 2.0 + (nlaps*sensorlength*sensorweight) / M3_norm_factor

    objective = M2 + M3



    ### --- Return outputs
    # we want to maximize objective, so need to swap sign upon return
    if return_all
        return -objective, M2, M3
    else
        return -objective
    end

end



function con2021!(con, x, p, c)

    ### --- Unpack Variables
    # Unpack applicable design variables
    cruisevelocity  = x[1]*10.0 # cruise velocity, feet/sec
    ncontainers     = x[2]*10.0 # number of containers
    # sensorlength    = x[3] # sensor length
    sensorweight    = x[4] # sensor weight
    batteryweight   = x[5] # battery weight
    batteryC        = x[6]*10.0 # battery C rating
    batteryvoltage  = x[7]*10.0 # battery voltage
    wingarea        = x[8] # wing area


    # Unpack applicable parameters
    laplength               = p[1] # assumed length of lap
    battery_specific_energy = p[4] # battery specific energy
    eta                     = p[5] # propulsive efficiency
    CLmax                   = p[6] # Wing CLmax
    rho                     = p[7] # Air Density
    emptyweightratio2       = p[8] # aircraft empty weight ratio M2
    emptyweightratio3       = p[9] # aircraft empty weight ratio M3
    gravity                 = p[10] # gravity acceleration
    LoverD                  = p[11] # Lift to Drag ratio
    CD0                     = p[12] # zero lift drag coefficient


    # Unpack applicable constraints
    ttotal              = c[1] # max time allowed for mission 3
    maxwingspan         = c[2] # max allowed wing span
    maxweight           = c[3] # max allowed total weight
    maxtakeoffdist      = c[4] # max takeoff distance
    maxbatterycapacity  = c[5] # max allowed battery power
    maxbendingmoment    = c[6]
    maxaspectratio      = c[7]



    ### --- Set up values to be constrained

    ## - Total Stored Battery Power
    battery_watthours = batteryweight*battery_specific_energy #watt-hours
    batterycapacity = battery_watthours/batteryvoltage #amp-hours
    grosspower = RCDesignSuite.battery_power(batterycapacity,batteryC,batteryvoltage) #watts



    ## - Weight
    # Get weight of everything besides battery and payload for each mission.
    # structuralweight2 = (batteryweight+sensorweight*ncontainers)*emptyweightratio2/(1-emptyweightratio2)
    # structuralweight3 = (batteryweight+sensorweight)*emptyweightratio3/(1-emptyweightratio3)
    wingweight = RCDesignSuite.wingweight(wingarea)

    # get total weight for missions 2 and 3
    weight2 = RCDesignSuite.sumweight([batteryweight;sensorweight*ncontainers;wingweight])
    weight3 = RCDesignSuite.sumweight([batteryweight;sensorweight;wingweight])



    ## - takeoff distance
    availablepower = RCDesignSuite.available_power(grosspower,eta)
    takeoffdist = RCDesignSuite.liftoffdistance(weight2,gravity,rho,wingarea,CLmax,availablepower*0.1)



    ## - Mission endurance
    endurance3 = RCDesignSuite.endurance_time(battery_specific_energy*3600, eta, LoverD, batteryweight, gravity, cruisevelocity, weight3)/60.0 * 1.25 #convert into minutes plus 25% safety factor



    ## - Double Checks
    # check endurance for mission 2 works
    tperlap = laplength/cruisevelocity
    # time for 3 laps on mission 2
    t2 = 3*tperlap * 2 # times 2 safety factor
    endurance2 = RCDesignSuite.endurance_time(battery_specific_energy*3600, eta, LoverD, batteryweight, gravity, cruisevelocity, weight2) # in seconds

    # check max velocity and cruise velocity work
    Ta = RCDesignSuite.available_thrust(availablepower, cruisevelocity)
    maxvel2 = RCDesignSuite.maxvelocity(Ta, weight2, wingarea, CD0)
    maxvel3 = RCDesignSuite.maxvelocity(Ta, weight3, wingarea, CD0)

    # Root bending moment:
    bendingmoment = RCDesignSuite.root_bending_moment(maxwingspan,weight2)

    ### --- Organize Constraints
    # stored power
    con[1] = (batterycapacity - maxbatterycapacity)/maxbatterycapacity

    # allowed weight
    con[2] = (weight2 - maxweight)/maxweight

    # takeoff distance
    con[3] = (takeoffdist - maxtakeoffdist)/maxtakeoffdist

    # sufficient endurance to last full time in mission 3
    con[4] = (ttotal - endurance3)/ttotal

    # sufficient endurance to accomplish 3 laps in mission 2
    con[5] = (t2 - endurance2)/t2

    # cruise velocity less than max for mission 2
    if maxvel2 == false
        con[6] = 1.0
    else
        con[6] = (cruisevelocity - maxvel2)/maxvel2
    end

    # cruise velocity less than max for mission 3
    if maxvel3 == false
        con[7] = 1.0
    else
        con[7] = (cruisevelocity - maxvel3)/maxvel3
    end


    # root bending moment can't break wing
    con[8] = (bendingmoment - maxbendingmoment)/maxbendingmoment

    # Aspect ratio isn't too big (replace with better root bending constraint above.)
    con[9] = (maxwingspan^2/wingarea - maxaspectratio)/maxaspectratio
end




function sensitivity2021(N=50)

    r = range(0.1,stop=1.9,length=N)

    x0, _, _, p, c = setup2021()

    obj = zeros(length(x0),N)
    M2 = zeros(length(x0),N)
    M3 = zeros(length(x0),N)

    maxm2 = 0.0
    maxm3 = 0.0

    for i=1:length(x0)

        desvar = copy(x0)

        for j=1:N

            desvar[i] = x0[i]*r[j]

            _, M2temp, M3temp = obj2021(desvar, p, c, return_all=true)

            if M2temp > maxm2
                maxm2 = M2temp
            end
            if M3temp > maxm3
                maxm3 = M3temp
            end

        end

    end

    p[2] = maxm2
    p[3] = maxm3

    for i=1:length(x0)

        desvar = copy(x0)

        for j=1:N

            desvar[i] = x0[i]*r[j]

            obj[i,j], M2[i,j], M3[i,j] = obj2021(desvar, p, c, return_all=true)

        end

    end

    obj0, M20, M30 = obj2021(x0, p, c, return_all=true)

    return -obj, M2, M3, -obj0, M20, M30, r

end

function intermediate_plots(r, obj, obj0;
                        labels = [
                            "Cruise Velocity";
                            "number of containers";
                            "sensor length";
                            "individual sensor weight";
                            "battery mass";
                            "battery C rating";
                            "battery voltage";
                            "wing area";])

    # Plot each variable on it's own plot so you can clearly see its sensitivity.
    for i=1:length(labels)
    figure()
    clf()
    plot(r,(obj[i,:].-obj0)./obj0,label=labels[i])
    legend()
    end

end

function final_plots(r, obj, obj0, labels;
    # TODO: identify the index of the variables you want
    toplot = [1;2;3;4],
    scalefactor = 1e3, # If needed, add a scale factor that helps the plot be more clear
    colors = ["C0";"C1";"C2";"C0";"C1";"C2"], # Pick your colors
    styles = ["-";"-";"-";"--";"--";"--"],  # Pick your styles.
    fs = (4.5,4.5*3/4), #figure size
    save_path = "./", #save path
    )

    r = r .- 1.0

    # Get the plot labels you want to show for those variables
    toplotlabels = labels[toplot]


    #### PLOT FINAL FIGURE FOR REPORTS ####
    # Initialize the figure so that it will be about the right size for the paper (4x3 ratio), but also such that the legend doesn't overlap things weird.
    figure(figsize=fs)

    xlabel("Percent Change in Design Variable (from nominal)")
    ylabel("Normalized Change in Score (x 1e$(Int(log10(scalefactor))))")

    # Plot the interesting things.
    for i=1:length(toplot)
    plot(r*1e2,scalefactor*(obj[toplot[i],:].-obj0)./obj0,label=toplotlabels[i],linestyle=styles[i],color=colors[i])
    end

    legend()

    #save the figure.
    savefig(save_path*"sensitivityobj.png",bbox_inches="tight")


end


obj, M2, M3, obj0, M20, M30, r = sensitivity2021(50)

intermediate_plots(r, obj, obj0)

labels = [
    "Cruise Velocity";
    "number of containers";
    "sensor length";
    "individual sensor weight";
    "battery mass";
    "battery C rating";
    "battery voltage";
    "wing area";]

final_plots(r, obj, obj0, labels)