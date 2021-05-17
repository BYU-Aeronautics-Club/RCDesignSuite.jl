
"""
Inputs
wing: the area of the wing
tail: the overall design of the tail
propulsion: the placement and design of the propulsion system
wieght: total weight of the aircraft
"""

struct concept
    wing            #Float64 describing the area of the wing in in^2
    tail            #String describing the tail configuration
    propulsion      #String describing placement and design for propulsion system
    weight          #Float64 describing the overal weight of the aircraft in lbs

end

aircraft = concept(350.0,"conventional","front propeller", 1.5)



"""
Preliminary aircraft type for storing prelminary design details
"""
struct prelim


end
