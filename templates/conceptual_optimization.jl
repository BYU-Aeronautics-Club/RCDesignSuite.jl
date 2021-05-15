#=

A first attempt at a conceptual  level optimization.

Authors: Judd Mehr,

=#


### --- Knowns --- ###

## -- Facts -- ##
g = 32.17405 # gravity, ft/s^2

## -- Constraints -- ##
M3_max_time = 10.0*60.0 # max time, in seconds, allowed for mission 3 flight.

wingspan = 5.0 # feet

maxweight = 55.0 # pounds

max_takeoff_distance = 100.0 # feet

max_turn_radius = 250.0 # feet, please don't turn this wide...

max_battery_energy = 200.0 # watt-hours






### --- Assumed Constants --- ###
battery_specific_energy = 35*2.2 # (watt-hours/kg) * (kg/pound) often called energy density, but this is a misnomer. #! Need to actually weight some typical batteries to estimate a good number here.

max_battery_weight = max_battery_energy/battery_specific_energy

propulsive_efficiency = 0.6 # guess. Never more than 0.8 for sure.

LoverD = 10.0 # lift/drag ratio. 10 = not excellent, but not absolutely terrible.


