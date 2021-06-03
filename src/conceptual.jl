#=

Authors: Judd Mehr, Andrew Tagg

Conceptual Design Functions

=#



############################
#######    General    ######
############################

"""
    endurance_time(eb, eta, LoverD, mb, g, Vinf, mto)

Calculate endurance time. (eqn. 7.21)

**Inputs:**

`eb::Float64` :

`eta::Float64` :

`LoverD::Float64` : Lift to Drag ratio, L/D.

`mb::Float64` :

`g::Float64` :

`Vinf::Float64` :

`mto::Float64` :
"""
function endurance_time(eb, eta, LoverD, mb, g, Vinf, mto)
    return (eb*eta*LoverD*mb)/(g*Vinf*mto)
end



"""
    climbrate(V,T,D,W)

Calculate climb rate. (eqn 7.52)

**Inputs:**

`V::Float64` : Freestream speed

`T::Float64` : Thrust

`D::Float64` : Drag

`W::Float64` : Weight
"""
function climbrate(V,T,D,W)
       return V*(T-D)/W
end




"""
    glideratio(L,D)

Calculate glide (lift to drag) ratio

**Inputs:**

`L::Float64` : Lift

`D::Float64` : Drag
"""
function glideratio(L,D)
    return L/D
end



"""
    turnradius(V,g,phi)

Calculate coordinated turn radius.

**Inputs:**
`V::Float64` : Freestream speed
`g::Float64` : Gravity
`phi::Float64` : Bank angle
"""
function turnradius(V,g,phi)
       return (V^2)/(g*tan(phi))
end



"""
    liftoffdistance(W,ρ,S,Cl,P)

Calculate distance from stand still to lift off. (eqn 7.81)

**Inputs:**

`W::Float64` : Weight

`ρ::Float64` : Denisty

`S::Float64` : Reference area

`Cl::Float64` : Maximum lift coefficient

`P::Float64` : Net power
"""
function liftoffdistance(W,ρ,S,Cl,P)
    return (1.629*W^(5/2))/(P*(ρ*S*Cl)^(3/2))
end



"""
    landingdistance(L,D,Vg,Vl,g)

Calculate deceleration distance before touchdown. (eqn 7.86)

**Inputs:**

`L::Float64` : Lift

`D::Float64` : Drag

`Vg::Float64` : Speed at the end of the glide phase

`Vl::Float64` : Speed immediately before touchdown

`g::Float64` : acceleration due to gravity
"""
function landingdistance(L,D,Vg,Vl,g)
    return (L/D)*((Vg^2)-(Vl^2))/(2*g)
end

"""
    landingdistance(Vl,g,mu,W,L,D)

Calculate deceleration distance after touchdown. (eqn 7.90)

**Inputs:**

`Vl::Float64` : Speed immediately before touchdown

`g::Float64` : acceleration due to gravity

`mu::Float64` : friction constant

`W::Float64` : Weight at landing

`L::Float64` : Lift

`D::Float64` : Drag
"""
function landingdistance(Vl,g,mu,W,L,D)
    return Vl^2/(2*g*(mu*(W-L)+D))
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
    vstall(weight,CLmax,rho,referencearea)

Calculate stall velocity.

**Inputs:**

`weight::Float64` : Aircraft weight.

`CLmax::Float64` : Aircraft maximum lift coefficient.

`rho::Float64` : Air density at flight conditions.

`referencearea::Float64` : Wing reference area.
"""
function vstall(weight,CLmax,rho,referencearea)
    return sqrt(weight/(CLmax*rho*referencearea))
end



"""
    vtaekoff(weight,CLmax,rho,referencearea)

Calculate takeoff velocity.

**Inputs:**

`weight::Float64` : Aircraft weight.

`CLmax::Float64` : Aircraft maximum lift coefficient.

`rho::Float64` : Air density at flight conditions.

`referencearea::Float64` : Wing reference area.
"""
function vtakeoff(weight,CLmax,rho,referencearea)
    return 1.2*sqrt(weight/(CLmax*rho*referencearea))
end

"""
    vtakeoff(vstall)

Calculate takeoff velocity.

**Inputs:**

`vstall::Float64` : Stall speed.
"""
function vtakeoff(vstall)
    return 1.2*vstall
end



"""
eqn 7.30
"""
function clmax_LD()

end



"""
    clmax_endurance(clmaxLD)

Calculate maximum endurance lift coefficient. (eqn 7.48)

**Inputs:**

`clmax_LD::Float64` : maximum lift to drag lift coefficient
"""
function clmax_endurance(clmaxLD)
    return sqrt(3)*clmaxLD
end




"""
    maxvelocity(Ta, W, S, rho, CD0, K=0.38)

Calculate maximum velocity.

**Inputs:**

`Ta::Float64` : maximum thrust available, in Newtons

`W::Float64` : total aircraft weight, in Newtons

`S::Float64` : area of wing in square meters

`CD0::Float64` : zero-lift drag coefficient of the aircraft.

`rho::Float64` : air density in kg/m^3, default is 1.225

`K::Float64` : empirical constant (default = 0.38)
"""
function maxvelocity(Ta, W, S, CD0, rho=1.225, K=0.38)

    num = ( (Ta/S) + (W/S)*sqrt((Ta/W)^2 - 4*CD0*K) )

    den = rho * CD0

    return sqrt( num / den )

end



"""
    oswald_factor(einv, CDp, AR, K=0.38)

Calculate the Oswald efficiency factor.

**Inputs:**

`einv::Float64` : inviscid span efficiency

`CDp::Float64` : parasitic drag coefficient

`AR::Float64` : wing aspectratio

`K::Float64` : empirical constant, default=0.38
"""
function oswald_factor(einv, CDp, AR, K=0.38)
    return 1 / ( (1/einv) + K*CDp*pi*AR)
end



"""
    eta_span(fuse_diameter, span)

Calculate inviscid span efficiency.

**Inputs:**

`fd::Float64` : fuselage diameter, in meters

`b::Float64` : wing span, in meters.
"""
function e(fd, b)
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

`capacity::Float64` : battery capacity in amp-hours

`C::Float64` : C-rating of LiPo battery

`voltage::Float64` : voltage of battery, in volts
"""
function battery_power(capacity, C, voltage)
    return capacity*C*voltage
end



"""
    available_power(gross_battery_power,eta_battery, eta_motor, eta_prop)

Calculate available power based on efficiencies.

**Inputs:**

`gross_battery_power::Float64` : theoretical max batter power, in Watts

`eta::Float64` : propulsive efficiency factor for battery, motor, and propeller (default is 0.8 which is the max theoretical, and high for real life)
"""
function available_power(gross_battery_power,eta=0.8)
    return gross_battery_power * eta
end



"""
    available_thrust(P, V)

Calculate thrust from available power and velocity.

**Inputs:**

`P::Float64` : power, in watts

`V::Float64` : velocity, in meters per seconds
"""
function available_thrust(P, V)
    return P/V
end


"""
    available_thrust(P, V)

Calculate thrust from available power and velocity.

**Inputs:**

`P::Float64` : power, in watts

`V::Float64` : velocity, in meters per seconds

`eta::Float64` : propulsive efficiency factor for battery, motor, and propeller
"""
function available_thrust(P, V, eta)
    return P*eta/V
end



############################
#####    Structures    #####
############################
"""
    turn_lift_force(weight,bankangle)

Calculates the lift force required in a coordinated turn.

**Inputs:**

`weight::Float64` : Weight of the aircraft.

`bankangle::Float64` : (radians) bank angle of aircraft.
"""
function turn_lift_force(weight,bankangle)
    return weight/cos(bankangle)
end



"""
    root_bending_moment(span,liftforce,liftcentroid)

Calculate root bending moment.

**Inputs:**

`span::Float64` : Total span of wing.

`liftforce::Float64` : Total  lift force on wing.

`liftcentroid::Float64` : point of equivalent point load. default = 4/(3*pi) for elliptic loading
"""
function root_bending_moment(span,liftforce,liftcentroid=4.0/(3.0*pi))
    return liftcentroid*span*liftforce/4.0
end



"""
    load_factor(lift,weight)

Calculate load factor, n.

**Inputs:**

`lift::Float64` : Lift Force

`weight::Float64` : Weight (in same units as lift)
"""
function load_factor(L,W)
    return L/W
end



"""
    equivalent_airspeed(trueairspeed, rho, rhosl)

Calculate equivalent air speed based on true airspeed, air density at flight condition, and sea level air density.

**Inputs:**

`trueairspeed::Float64` : free stream velocity at flight condition.

`rho::Float64` : air density at flight condition.

`rhosl::Float64` : air density at sea level. default = 1.225 kg/m^3
"""
function equivalent_airspeed(trueairspeed, rho, rhosl=1.225)
    return trueairspeed*sqrt(rho/rhosl)
end



"""
    Vn_stall(equivalentairspeed, referencearea, weight, CLmax, rhosl)

Get load factor value along stall curve for V-n diagram as a function of equivalent airspeed.

**Inputs:**

`equivalentairspeed::Float64` : Equivalent Air Speed (velocity, i.e. the x-axis value of V-n diagram)

`referencearea::Float64` : Wing reference area.

`weight::Float64` : Total weight of aircraft.

`CLmax::Float64` : Maximum Lift Coeficient of the aircraft. default = 1.2

`rhosl::Float64` : air density at sea level. default = 1.225 kg/m^3
"""
function Vn_stall(equivalentairspeed, referencearea, weight, CLmax=1.2, rhosl=1.225)
    return CLmax*rhosl*equivalentairspeed^2*referencearea/(2*weight)
end



"""
    weight(weights)

Sums up weights and returns total.

**Inputs:**

`weight::Array{Float64}` : Array of weights
"""
function weight(weights)
    return sum(weights)
end
