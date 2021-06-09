"""
Conceptual aircraft type for storing conceptual design details
"""
struct concept


end



"""
Preliminary aircraft type for storing prelminary design details
"""
struct prelim
    struct prelim
        span
        thickness
        chord
        twist angle
        dihedral
        angle of attack
        sideslip angle
        bank angle
        ns                  #number of panels along span, for use with vortex lattice
        nc                  #number of panels along chord
    end

end
