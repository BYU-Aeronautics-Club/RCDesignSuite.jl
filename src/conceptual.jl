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

`eb::Float64` : battery specific energy

`eta::Float64` : propulsive efficiency factor

`LoverD::Float64` : Lift to Drag ratio, L/D.

`mb::Float64` : battery mass

`g::Float64` : gravity

`Vinf::Float64` : cruise velocity

`mto::Float64` : total aircaft mass at takeoff
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

`rho::Float64` : Air density

`S::Float64` : Reference area

`Clmax::Float64` : Maximum lift coefficient

`P::Float64` : Net power including drag
"""
function liftoffdistance(W,g,rho,S,Clmax,P)
    return (1.629*W^(5/2))/(g*P*(rho*S*Clmax)^(3/2))
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


"""
    dragparasitic(S, v; tc = .12, AR = 6, fr = 4, Λ = 4, ρ = 1.225, μ = 1.81e-5)

Calculates the sum of the parasitic drag on the fuselage and on the wing

**Inputs:**

`S::Float64` : reference area

`v::Float64` : velocity

`tc::Float64` : wing thickness to chord ratio - thickness to chord ratio assumed to be about .12

`AR::Float64` : aspect ratio assumed to be about 6

`fr::Float64` : fineness ratio assumed to be about 4

`Λ::Float64` : sweep angle assumed to be about 4 degrees

`ρ::Float64` : air density assumed to be 1.225 kg/m^3

`μ::Float64` : air viscosity assumed to be 1.81e-5 kg/ms
"""
function dragparasitic(S, v; tc = .12, AR = 6, fr = 4, Λ = 4, ρ = 1.225, μ = 1.81e-5)
    q = .5*ρ*v^2
    b = sqrt(AR*S)
    c = S/b
    return dragwing(q, S, v, c, Λ, tc, ρ, μ) + dragfuselage(q, S, v, c, fr, ρ, μ)
end


"""
    dragwing(S, v; tc = .12, AR = 6, fr = 4, Λ = 4, ρ = 1.225, μ = 1.81e-5)

Calculates the parasitic drag on the wing as part of the dragparasitic function

**Inputs:**

`q::Float64` : dynamic pressure

`S::Float64` : reference area

`v::Float64` : velocity

`c::Float64` : chord

`Λ::Float64` : sweep angle

`tc::Float64` : wing thickness to chord ratio - thickness to chord ratio

`ρ::Float64` : air density

`μ::Float64` : air viscosity
"""
function dragwing(q, S, v, c, Λ, tc, ρ, μ)
    z = 2*cos(Λ)
    k = 1 + z*(tc) + 100*(tc)^4
    Re = ρ*v*c/μ
    Cf = .074/(Re^.2)
    Swet = 2*(1+.2*(tc))*S
    Dp = k*Cf*q*Swet
    return Dp
end


"""
    dragfuselage(S, v; tc = .12, AR = 6, fr = 4, Λ = 4, ρ = 1.225, μ = 1.81e-5)

Calculates the parasitic drag on the fuselage as part of the dragparasitic function

**Inputs:**

`q::Float64` : dynamic pressure

`S::Float64` : reference area

`v::Float64` : velocity

`c::Float64` : chord

`fr::Float64` : fineness ratio assumed to be about 4

`ρ::Float64` : air density

`μ::Float64` : air viscosity
"""
function dragfuselage(q, S, v, c, fr, ρ, μ)
    l = S/c
    k = 1.675 - 0.09*fr+0.003*fr^2
    if fr >= 15
        k = 1
    end
    Re = ρ*v*l/μ
    Cf = .074/(Re^.2)
    S = pi*(l^2/fr)*1.75
    Dp = k*Cf*q*S
    return Dp
end


"""
    draginduced(w, S, v; ρ = 1.225)

Calculates the induced drag needed to produce lift

**Inputs**

`w::Float64` : weight

`S::Float64` : reference area

`v::Float64` : velocity

`ρ::Float64` : air density assumed to be 1.225 kg/m^3
"""
function draginduced(w, S, v; ρ = 1.225)
    q = .5*ρ*v^2
    CL = w/(q*S)
    Dp = dragparasitic(S, v)
    Cp = Dp/(q*S)
    Ci = .38*Cp*(CL^2)
    Di = Ci*q*S
    return Di
end


"""
    dragdata(a::concept)

Outputs the induced, parasitic, and total drag of an aircraft with respect to velocity

**Inputs**

A struct of type `concept`
"""
function dragdata(a::concept)
    w = a.weight
    S = a.area
    dragtot = []
    dragpar = []
    dragin =[]
    velocity = []
    vs = range(1, stop = 50, length = 201)
    for v in vs
        Dp = dragparasitic(S, v)
        Di = draginduced(w, S, v)
        Dt = Dp+Di
        push!(dragtot, Dt)
        push!(dragpar, Dp)
        push!(dragin, Di)
        push!(velocity, v)
        v += 1
    end

    return velocity, dragpar, dragin, dragtot
end


"""
    dragcalculator(a::concept, v)

Outputs the total drag on an aircraft

**Inputs**

A struct of type `concept`

`v::Float64` : velocity
"""
function dragcalculator(a::concept, v)
    w = a.weight
    S = a.area
    Dp = dragparasitic(S, v)
    Di = draginduced(w, S, v)
    Dt = Dp+Di
    return Dt
end


"""
    designspeed(a::concept)
Calculates the most efficient speed of the aircraft based on drag

**Inputs**

A struct of type `concept`
"""
function designspeed(a::concept)
    w = a.weight
    S = a.area
    vs = range(1, stop = 50, length = 201)
    Dt = zeros(length(vs))
    for (i,v) in enumerate(vs)
        Dp = dragparasitic(S, v)
        Di = draginduced(w, S, v)
        Dt[i] = Dp + Di
    end
    minimumdrag,iminimumdrag = findmin(Dt)
    return vs[iminimumdrag]

end


"""
    liftcoefficient(a::concept;ρ = 1.225)
Calculates the data for lift coefficient with respect to velocity.  Outputs the design lift coefficient based on the design speed

**Inputs**

A struct of type `concept`
"""
function liftcoefficient(a::concept;ρ = 1.225)
    w = a.weight
    S = a.area
    vs = range(1, stop = 50, length = 201)
    CL = zeros(length(vs))
    designliftcoeff= 0.0
    vdesign = designspeed(a)
    for (i, v) in enumerate(vs)
        q = .5*ρ*v^2
        CL[i] = w/(q*S)
        if v == vdesign
            designliftcoeff = CL[i]
        end
    end
    println("design lift coefficient is: ", designliftcoeff)
    return vs, CL, designliftcoeff
end



######################################################################
############################    Sizing    ############################
######################################################################
"""
"""
function required_wing_area()

end


"""
    get_fuselage_volume(n,v)

Calculate fuselage volume based on number and volume of repeated payload components.  Also optional argument to add miscellanous extra volume items (batteries, ESC's, etc.)

**Inputs:**

`n::Int` : number of repeated payload items (internally stored in fuselage)

`v::Float64` : volume of singe payload item.

`miscvol::Float64` : additional volume, e.g. from batteries, esc's, etc.
"""
function get_fuselage_volume(n,v,miscvol=0.0)
    return v*n + miscvol
end


######################################################################
#########################    Aerodynamics    #########################
######################################################################

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
    clmax_LD(CDp,AR,e)

Calculate lift coefficient for maximum L/D (eqn 7.30)

**Inputs:**

`CDp::Float64` : Parasitic drag coefficient

`AR::Float64` : Aspect ratio

`e::Float64` : Oswald efficiency factor
"""
function clmax_LD(CDp,AR,e)
    return sqrt(CDp*pi*AR*e)
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

    if ((Ta/W)^2 - 4*CD0*K) < 0.0
        return false
    else
        num = ( (Ta/S) + (W/S)*sqrt((Ta/W)^2 - 4*CD0*K) )

        den = rho * CD0

        return sqrt( num / den )
    end

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



#######################################################################
###########################    Propulsion    ##########################
#######################################################################

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



########################################################################
###########################    Structures    ###########################
########################################################################
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
function sumweight(weights)
    return sum(weights)
end


"""
    wingweight()

Estimates Wing Weight.

**Inputs:**

`wingarea::Float64` : wing area.

**Keyword Arguments:**
`volume::Float64` : volume

`density::Float64` : density, , default is from lookup table for pink foam in kg/m^3

`foamdatafile::String` : if using, path to foam data file, default = "./data/foamdata_SI.csv",

`foamtype::String` : if using, FoamType to use in lookup table, default = "Foamular250a" (pink foam)

`afthickness::Float64` : airfoil percent thickness, default = 0.14

`span::Float64` : wing span, default = 1.524 m
"""
function wingweight(wingarea;
            volume          = nothing,
            density         = nothing,
            foamdatafile    = "./data/foamdata_SI.csv",
            foamtype        = "Foamular250a",
            afthickness     = 0.14,
            span            = 1.524)

    # if density isn't manually input, get it from lookup table.
    if density == nothing

        foamdata = DataFrames.DataFrame(CSV.File(foamdatafile,header=true,datarow=3))

        foamidx = findall((foamdata.FoamType .== foamtype))

        density = foamdata.rho[foamidx] #kg/m^3

    end

    # if volume isn't manuual input, estimate it based on span and thickness.
    if volume == nothing

        # get average chord
        cbar = wingarea/span

        #multiply current thickness by the ratio of area to thickness for NACA 2412
        afarea = cbar*afthickness*(0.081699375/0.12)

        # multiply airfoil area by wing area
        wingvolume = afarea*span

    end

    # Return weight
    return wingvolume*density

end



"""
Need Method for coupling wing area and root bending yield stress.
Probabaly need to have a parameter or variable for airfoil thickness
Get chord from area and span parameter
Get actual thickness from chord and thickness parameter/variable
Calculate bending yield stress as a function of absolute section thickness.
"""
function bending_yield_stress()

end