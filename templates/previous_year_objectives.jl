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
        50.0    # Cruise Velocity, feet/second
        10.0    # number of containers
        1.0     # sensor length, feet
        0.125   # individual sensor weight, pounds
        0.5     # battery mass, pounds
        25.0    # battery C rating
        21.0    # battery voltage
        2.5     # wing area
    ]

    # Lower Bounds
    lb = [
        0.0     # Cruise Velocity, feet/second
        1.0     # number of containers
        0.01    # sensor length, feet
        0.01    # individual sensor weight, pounds
        0.0     # battery mass, pounds
        10.0    # battery C rating
        3.7     # battery voltage
        0.0     # wing area
    ]

    # Upper Bounds
    ub = [
        150.0   # Cruise Velocity, feet/second
        100.0   # number of containers
        2.0     # sensor length, feet
        5.0     # individual sensor weight, pounds
        2.0     # battery mass, pounds
        50.0    # battery C rating
        34.0    # battery voltage
        25.0    # wing area
    ]

    # Parameter Values
    p = [
        4000.0  # assumed length of lap, feet
        1.0     # normalization factor for mission 2
        1.0     # normalization factor for mission 3
    ]

    # Constraint Values
    c = [
        600.0   # max time, in seconds, allowed for mission 3 flight.
        5.0     # max wingspan, feet (projected)
        55.0    # maximum allowed weight, pounds
        100.0   # maximum takeoff distance, feet
        200.0   # maximum battery energy watt-hours
    ]

    return x0, lb, ub, p, c

end



function obj2021(x,p, c)

    ### --- Unpack Variables
    # Unpack applicable design variables
    cruisevelocity  = x[1] # cruise velocity, feet/sec
    ncontainers     = x[2] # number of containers
    sensorlength    = x[3] # sensor length
    sensorweight    = x[4] # sensor weight

    # Unpack applicable parameters
    laplength       = p[1] # assumed length of lap
    M2_norm_factor  = p[3] # normalization factor for mission 2
    M3_norm_factor  = p[4] # normalization factor for mission 3

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



function con2021()

end