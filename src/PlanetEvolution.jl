module PlanetEvolution

greet() = print("Hello Explorer!")

include("tidal.jl")
using .TRIPS

include("thermal.jl")
using .ThE
export
    one_layer_plnt,
    two_layer_plnt,
    find_core

include("orbits.jl")
using .Orbit
export
    orbital_recession


end # module
