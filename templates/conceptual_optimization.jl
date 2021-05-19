#=

A first attempt at a conceptual  level optimization.

Authors: Judd Mehr,

=#

#############################################
##########          SETUP          ##########
#############################################

### --- KNOWNS --- ###

## -- Facts -- ##
g = 32.17405 # gravity, ft/s^2
rhosl = 0.0765 # air density, sea level  lb/ft^3 = 1.225 kg/m^3

## -- Constraints -- ##
M3_max_time = 10.0*60.0 # max time, in seconds, allowed for mission 3 flight.

wingspan = 5.0 # feet (projected)

maxweight = 55.0 # pounds

max_takeoff_distance = 100.0 # feet

max_battery_energy = 200.0 # watt-hours






### --- ASSUMED CONSTANTS (Parameters) --- ###
battery_specific_energy = 35*2.2 # (watt-hours/kg) * (kg/pound) often called energy density, but this is a misnomer. #! Need to actually weigh some typical batteries to estimate a good number here.

max_battery_weight = max_battery_energy/battery_specific_energy

propulsive_efficiency = 0.6 # guess. Never more than 0.8 for sure.

LoverD = 10.0 # lift/drag ratio. 10 = not excellent, but not absolutely terrible.

rho = rhosl # assume sea level flight until picking a flight test location.

aircraft_empty_weight = 0.0 # pounds, need an experienced-based number, perhaps as a function of total weight.

CD0 = 0.0 # zero-lift drag coefficient of aircraft, need educated guess.

lap_length = 2*500*(1+pi) # somewhat arbitrary, really big turns probably account for extra loop length.




### --- VARIABLES --- ###
battery_mass = 0.5 # pounds

battery_C_rating = 25.0 # need this to get available power

battery_voltage = 21.0 # need this to get available power

wing_area = 2.5 # square feet

payload_weight = 0.0 # pounds, need to know payload specific requirements

Vinf = 75.0 # cruise velocity, feet/sec



#############################################
##########        OBJECTIVE        ##########
#############################################

function obj()

end