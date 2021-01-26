using PyPlot
using Roots
using FLOWMath


##########################################################################
####################         MISSION OBJECTIVES       ####################
##########################################################################

"""
    mission2(inputs)

Calculates score from mission 2 based on relavent inputs.

**Inputs:**
- num_sensors::Int64 : the number of sensors.
- time::Float64 : the time to complete 3 laps, in seconds.
"""
function mission2(num_sensors, time)

    global maxM2 #note, using globals like this is poor form, but convenient.

    M2 = num_sensors/time

    if M2 > maxM2
        global maxM2 = M2
    end

    return 1.0 + M2/maxM2
end


"""
    mission3(inputs)

Calculates score from mission 3 based on relavent inputs.

**Inputs:**
- num_laps::Int64 : the number of laps completed in 10 minutes.
- length_sensor::Float64 : the length of the sensor, in meters.
- weight_sensor::Float64 : the weight of the sensor, in kilograms.
"""
function mission3(num_laps, length_sensor, weight_sensor)
    global maxM3

    M3 = num_laps * length_sensor * weight_sensor


    if M3 > maxM3
        global maxM3 = M3
    end

    return 2.0 + M3/maxM3
end



##########################################################################
####################         HAND CALCULATIONS        ####################
##########################################################################

"""
    get_thrust(P, V)

Calculate thrust from available power and velocity.

**Inputs:**
- P::Float64 : power, in watts
- V::Float64 : velocity, in meters per seconds
"""
function get_thrust(P, V)
    return P/V
end


"""
    gross_battery_power(capacity, C, voltage)

Calculate theoretical max battery power.

**Inputs:**
- capacity::Float64 : battery capacity in amp-hours
- C::Float64 : C-rating of LiPo battery
- voltage::Float64 : voltage of battery, in volts
"""
function get_battery_power(capacity, C, voltage)
    return capacity*C*voltage
end


"""
    available_power(gross_battery_power,eta_battery, eta_motor, eta_prop)

Calculate available power based on efficiencies.

**Inputs:**
- gross_battery_power::Float64 : theoretical max batter power, in Watts
- eta::Float64 : efficiency factor for battery, motor, and propeller (propulsive efficiency factor)
"""
function get_available_power(gross_battery_power,eta)
    return gross_battery_power * eta
end


"""
    maxvelocity(Ta, W, S, rho, CD0, K=0.38)

Calculate maximum velocity.

**Inputs:**
- Ta::Float64 : maximum thrust available, in Newtons
- W::Float64 : total aircraft weight, in Newtons
- S::Float64 : area of wing in square meters
- CD0::Float64 : zero-lift drag coefficient of the aircraft.
- rho::Float64 : air density in kg/m^3, default is 1.225
- K::Float64 : empirical constant (default = 0.38)
"""
function maxvelocity(Ta, W, S, CD0, rho=1.225, K=0.38)

    num = ( (Ta/S) + (W/S)*sqrt((Ta/W)^2 - 4*CD0*K) )

    den = rho * CD0

    return sqrt( num / den )

end


"""
    get_weight(weights)

Sums up weights and returns total.

**Inputs:**
- weight::Array{Float64} : Array of weights
"""
function get_weight(weights)
    return sum(weights)
end

#= Currently Unused
"""
    get_oswald_factor(einv, CDp, AR, K=0.38)

Calculate the Oswald efficiency factor.

**Inputs:**
- einv::Float64 : inviscid span efficiency
- CDp::Float64 : parasitic drag coefficient
- AR::Float64 : wing aspectratio
- K::Float64 : empirical constant, default=0.38
"""
function get_oswald_factor(einv, CDp, AR, K=0.38)
    return 1 / ( (1/einv) + K*CDp*pi*AR)
end


"""
eta_span(fuse_diameter, span)

Calculate inviscid span efficiency.

**Inputs:**
- fd::Float64 : fuselage diameter, in meters
- b::Float64 : wing span, in meters.
"""
function get_e(fd, b)
    return 0.98*(1-2*(fd/b)^2)
end
=#

############################################################################
#####################        OBJECTIVE FUNCTION        #####################
############################################################################

function objective(x0)

    #unpack design parameters
    W               = x0[1]
    S               = x0[2]
    CD0             = x0[3]
    lap_distance    = x0[4]
    num_sensors     = x0[5]
    length_sensor   = x0[6]
    weight_sensor   = x0[7]
    M2_laps         = x0[8]
    maxtime         = x0[9]
    eta             = x0[10]
    batteryC        = x0[11]
    batterycells    = x0[12]
    batterycapacity = x0[13]

    batteryvoltage = batterycells*3.7

    P = get_battery_power(batterycapacity, batteryC, batteryvoltage)
    Pa = get_available_power(P,eta)

    ### MISSION 2 ###

    #get loaded weight for mission 2
    W2 = get_weight([W; weight_sensor*num_sensors])

    #get velocity for mission 2, assume max velocity the whole time
    findV2(V) = maxvelocity(get_thrust(Pa,V), W2, S, CD0) - V
    V2 = fzero(findV2,1.0,1000.0)[1]

    #get time for mission 2 completion
    time2 = lap_distance*M2_laps/V2


    ### MISSION 3 ###

    #get loaded weight for mission 3
    W3 = get_weight([W; weight_sensor])

    #get velocity for mission 3, assume max velocity the whole time
    findV3(V) = maxvelocity(Pa/V, W3, S, CD0) - V
    V3 = find_zeros(findV3,1.0,1000.0)[1]

    #get time for a single lap in mission 3
    time3 = lap_distance/V3

    #get number of laps completable in time
    num_laps = maxtime/time3

    ### SUM UP OBJECTIVES ###

    M2 = mission2(num_sensors, time2)

    #convert units for metrics:
    length_sensor *= 39.3701 #meters to inches
    weight_sensor *= 35.274 #kg to oz

    M3 = mission3(num_laps, length_sensor, weight_sensor)
    totalpoints = 100.0 + 1.0 + 1.0 + M2 + M3 #assumes perfect on report and GM and M1 completion.

    return totalpoints

end



############################################################################
#####################        SENSITIVITY STUDY         #####################
############################################################################

#for missions 2 and 3, create globals so we can find the max of our sensitivity study in order to have a denominator for the overall scoring, keeping sensitivities in the same order of magnitude.
global maxM2 = 0.0
global maxM3 = 0.0


#Define all the variables that might be important.
W               = 2.5       # aircraft weight (including battery), kilograms
S               = 2.5       # aircraft wing area, meters squared
CD0             = 0.01      # zero-lift parasitic drag coefficient
lap_distance    = 1000.0    # lap distance, meters, from rules image: (500.0*4.0 + pi*285.0 + pi*200.0)*0.3048
num_sensors     = 5.0       # number of sensors
length_sensor   = 0.15      # length of sensors, meters
weight_sensor   = 0.1       # weight of sensors, kilograms
M2_laps         = 3         # laps for mission 2
maxtime         = 10*60     # maximum time for mission 3
eta             = 0.75      # efficiency of propulsion system (battery, motor, prop)
batteryC        = 25        # C-rating of battery
batterycells    = 4         # number of battery cells (assumes single battery)
batterycapacity = 3         # battery capacity (assumes single battery)


#Define input array for objective function.
x0      = zeros(13)
x0[1]   = W
x0[2]   = S
x0[3]   = CD0
x0[4]   = lap_distance
x0[5]   = num_sensors
x0[6]   = length_sensor
x0[7]   = weight_sensor
x0[8]   = M2_laps
x0[9]   = maxtime
x0[10]  = eta
x0[11]  = batteryC
x0[12]  = batterycells
x0[13]  = batterycapacity

#define range of percentages to do sensitivity study over.
r = range(-50,stop=50,length=50)

#number of steps in the range
N = 50 #number of points in ranges

#initialize the objectives.
obj = zeros(length(x0),N)
#initialize the objective derivatives.
dobj = zeros(length(x0),N)

#for each input, loop through a +/- 50% range from the nominal
#Run once to get maxima for maxM2 and maxM3 globals for denominators.
for i = 1:length(x0)
    x0copy = copy(x0)
    for j = 1:N
        x0copy[i] = (1.0 + 1e-2*r[j])*x0[i]
        #try catch in case you get a combination that doesn't make sense or has a negative square root or something.
        try
            #throw away the answer, we don't need it, we're just getting the globals to be the right values.
            _ = (objective(x0copy)-obj0)/obj0
        catch
            println("i: $i, j: $j")
        end
    end
end

#run again to get actual values for the objectives, normalized by the max M2 and M3 values just found.
for i = 1:length(x0)

    x0copy = copy(x0)

    for j = 1:N

        x0copy[i] = (1.0 + 1e-2*r[j])*x0[i]
        try
            obj[i,j] = objective(x0copy)
        catch
            obj[i,j] = NaN
        end
    end

    #spline the outputs.
    so = Akima(r, obj[i,:])

    #get the derivatives relative to the percent of input.
    for j = 1:N
        dobj[i,j] = derivative(so,r[j])
    end

end

#get the objective value for the nominal case for getting the percent differences later.
obj0 = objective(x0)


############################################################################
#####################        PLOTS TO LOOK AT          #####################
############################################################################

#create a vector for the labels, so you know what the following plots are for.
labels = [
"W";
"S";
"CD0";
"lap distance";
"\\# sensors";
"sensor length";
"sensor weight";
"M2 laps";
"maxtime";
"eta";
"batteryC";
"batterycells";
"batterycapacity"]

#plot each variable on it's own plot so you can clearly see its sensitivity.
for i=1:length(x0)
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

#identify the index of the variables you want
#=
Note: in this example, all the payload things have the same sensitivities, so we'll put them all together.  Also the wing loading and drag coefficient are the same, so they are together.  It also turns out that weight doesn't matter, so we'll  leave it out.  Furthermore, all the battery and efficiency stuff goes into available power, so we'll combine them under that heading.
=#
toplot = [1;2;3;5;10]
#get the plot labels you want to show for those variables
toplotlabels = ["Weight";"Wing Area"; "Parasitic Drag"; "Payload"; "Available Power"]
#pick your colors
colors = ["C0";"C1";"C3";"C0";"C1";"C2"]
#pick your styles.
styles = ["-";"-";"-";"--";"--";"--"]



#### PLOT FINAL FIGURE FOR PAPER ####
#initialize the figure so that it will be about the right size for the paper (4x3 ratio), but also such that the legend doesn't overlap things weird.
figure(figsize=(4.5,4.5*3/4))
clf()
xlabel("Percent Change in Design Variable (from nominal)")
ylabel("Normalized Change in Score")
#plot the interesting things.
for i=1:length(toplot)
    plot(r,1e3*(obj[toplot[i],:].-obj0)./obj0,label=toplotlabels[i],linestyle=styles[i],color=colors[i])
end
legend()
#save the figure.
savefig("../reports/template/designreport/figures/sensitivityobj.pdf",bbox_inches="tight")


#### MAYBE PLOT THE DERIVATIVES INSTEAD, IF THAT MAKES MORE SENSE ####
#initialize the figure so that it will be about the right size for the paper (4x3 ratio), but also such that the legend doesn't overlap things weird.
figure(figsize=(4.5,4.5*3/4))
clf()
xlabel("Percent Change in Design Variable (from nominal)")
ylabel("Differential Percent Change in Score (dScore/dx)")
#plot the interesting things.
for i=1:length(toplot)
    plot(r,dobj[toplot[i],:],label=toplotlabels[i],linestyle=styles[i],color=colors[i])
end
legend()
#save the figure.
savefig("../reports/template/designreport/figures/sensitivitydobj.pdf",bbox_inches="tight")
