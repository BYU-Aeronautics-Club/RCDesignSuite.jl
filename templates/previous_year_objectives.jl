#=

Previous Years setup, objective, and constraint functions for examples and comparisons.

Authors: Judd Mehr,

=#


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
        0.0     # battery mass, pounds
        10.0/10.0    # battery C rating
        3.7/10.0     # battery voltage
        0.0     # wing area , m^2
    ]

    # Upper Bounds
    ub = [
        50.0/10.0   # Cruise Velocity, m/s
        100.0/10.0   # number of containers
        1.0     # sensor length, meters
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
        100.0   # battery specific energy, watt-hours/kg
        0.6     # propulsive efficiency
        0.8     # wing CL max
        1.225   # air density, kg/m^3
        0.50    # M2 empty weight ratio (ratio of empty weight to total weight)
        0.75    # M3 empty weight ratio
        9.81    # gravity acceleration constant
        10.0    # lift to drag ratio
    ]

    # Constraint Values
    c = [
        10.0   # max time, in minutes, allowed for mission 3 flight.
        1.524   # max wingspan, meters (projected)
        25.0    # maximum allowed weight, kg
        30.0    # maximum takeoff distance, meters
        200.0   # maximum battery stored power watt-hours
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
    # laplength               = p[1] # assumed length of lap
    battery_specific_energy = p[4] # battery specific energy
    eta                     = p[5] # propulsive efficiency
    CLmax                   = p[6] # Wing CLmax
    rho                     = p[7] # Air Density
    emptyweightratio2       = p[8] # aircraft empty weight ratio M2
    emptyweightratio3       = p[9] # aircraft empty weight ratio M3
    gravity                 = p[10] # gravity acceleration
    LoverD                  = p[11] # Lift to Drag ratio


    # Unpack applicable constraints
    ttotal          = c[1] # max time allowed for mission 3
    # maxwingspan     = c[2] # max allowed wing span
    maxweight       = c[3] # max allowed total weight
    maxtakeoffdist  = c[4] # max takeoff distance
    maxbatterycapacity = c[5] # max allowed battery power



    ### --- Set up values to be constrained

    ## - Total Stored Battery Power
    battery_watthours = batteryweight*battery_specific_energy #watt-hours
    batterycapacity = battery_watthours/batteryvoltage #amp-hours
    grosspower = battery_power(batterycapacity,batteryC,batteryvoltage) #watts



    ## - Weight
    # Get weight of everything besides battery and payload for each mission.
    structuralweight2 = (batteryweight+sensorweight*ncontainers)*emptyweightratio2/(1-emptyweightratio2)
    structuralweight3 = (batteryweight+sensorweight)*emptyweightratio3/(1-emptyweightratio3)

    # get total weight for missions 2 and 3
    weight2 = weight([batteryweight;sensorweight*ncontainers;structuralweight2])
    weight3 = weight([batteryweight;sensorweight;structuralweight3])



    ## - takeoff distance
    availablepower = available_power(grosspower,eta)
    takeoffdist = liftoffdistance(weight2,gravity,rho,wingarea,CLmax,availablepower)



    ## - Mission 3 endurance
    endurance = endurance_time(battery_specific_energy, eta, LoverD, batteryweight, gravity, cruisevelocity, weight3)*60.0 #convert into minutes



    ### --- Organize Constraints
    con[1] = (batterycapacity - maxbatterycapacity)/maxbatterycapacity # stored power
    con[2] = (weight2 - maxweight)/maxweight # allowed weight
    con[3] = (takeoffdist - maxtakeoffdist)/maxtakeoffdist # takeoff distance
    con[4] = (ttotal - endurance)/ttotal # sufficient endurance
end