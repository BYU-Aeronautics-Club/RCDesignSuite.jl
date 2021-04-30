#=

Authors: Judd Mehr,

Conceptual Design functions

=#



############################
#######    General    ######
############################

"""
eqn 7.21 in 415 book
"""
function endurancetime()

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

############################
#####    Propulsion    #####
############################

#! See section 7.1.2,  in 415 book


############################
#####    Structures    #####
############################

#! See chapter 8 in book