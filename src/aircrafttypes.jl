
"""
Conceptual aircraft type for storing conceptual design details

Inputs

* `wing::Float64`: describes the area of the wing in m^2
* `tail::String`: describes the tail configuration
* `propulsion::String`: describes placement and design for propulsion system
* `weight::Float64`: describes the overall weight of the aircraft in kg
"""
struct concept
    wing::Float64           #describes the area of the wing in m^2
    tail::String            #describes the tail configuration
    propulsion::String      #describes placement and design for propulsion system
    weight::Float64         #describes the overall weight of the aircraft in kg

end


"""
Preliminary aircraft type for storing prelminary design details
"""
struct prelim


end
