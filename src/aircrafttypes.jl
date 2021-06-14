"""
Conceptual aircraft type for storing conceptual design details
"""
struct concept{TF, TS}
    area::TF           #describes the area of the wing
    weight::TF         #total weight of the aircraft
    tail::TS            #describes the design of the tail i.e. "conventional"
    propulsion::TS      #describes the propulsion system of the aircraft
end



"""
Preliminary aircraft type for storing prelminary design details
"""
    struct prelim{TF}
        span::TF
        thickness::TF
        rootchord::TF
        tipchord::TF
        twistangle::TF
        dihedral::TF
        polyhedral::TF
        angleofattack::TF
        sweep::TF
        aspectratio::TF
        taperratio::TF
        finenessratio::TF
        cg::TF
    end
