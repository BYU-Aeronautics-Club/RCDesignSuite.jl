#=

Authors: Judd Mehr,

Conceptual Design functions

=#



############################
#######    General    ######
############################

"""
eqn 7.21 in 415 book (need to fill out docstring)
"""
function endurance_time(eb, eta, L, mb, g, Vinf, D, mto)
    return (eb*eta*L*mb)/(g*Vinf*D*mto)
end



"""
eqn 7.52
"""
function climbrate()

end



"""
eqn 7.57
"""
function glideratio()

end



"""
eqn 7.60
"""
function turnradius()

end



"""
eqn 7.82 (also note potential need for climb phase models)
"""
function liftoffdistance()

end



"""
eqns 7.87 and 7.91?
"""
function landingdistance()

end



############################
#######    Sizing    #######
############################
"""
"""
function wingareareq()

end


############################
####    Aerodynamics    ####
############################

"""
eqn 4.1 in 415 book
"""
function vstall()

end



"""
"""
function vtakeoff(vstall)
    return 1.2*vstall
end



"""
eqn 7.30
"""
function clopt()

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



############################
#####    Propulsion    #####
############################

#! See section 7.1.2,  in 415 book for relevant functions



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
- eta::Float64 : propulsive efficiency factor for battery, motor, and propeller (default is 0.8 which is the max theoretical, and high for real life)
"""
function get_available_power(gross_battery_power,eta=0.8)
    return gross_battery_power * eta
end



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



############################
#####    Structures    #####
############################

#! See chapter 8 in book



"""
    get_weight(weights)

Sums up weights and returns total.

**Inputs:**
- weight::Array{Float64} : Array of weights
"""
function get_weight(weights)
    return sum(weights)
end