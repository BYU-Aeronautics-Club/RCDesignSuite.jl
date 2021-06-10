"""
Conceptual aircraft type for storing conceptual design details
"""
struct concept
    area::Float64           #describes the area of the wing
    weight::Float64         #total weight of the aircraft
    tail::String            #describes the design of the tail i.e. "conventional"
    propulsion::String      #describes the propulsion system of the aircraft
end



"""
Preliminary aircraft type for storing prelminary design details
"""
    struct prelim
        span::Float64
        thickness::Float64
        rootchord::Float64
        tipchord::Float64
        twistangle::Float64
        dihedral::Float64
        polyhedral::Float64
        angleofattack::Float64
        sweep::Float64
        aspectratio::Int64
        taperratio::Float64
        finenessratio::Float64
        cg::Float64

    end
