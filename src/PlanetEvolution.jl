module PlanetEvolution

greet() = print("Hello Explorer!")

include("util.jl")
using .Util
export
    write_json,
    read_json

include("tidal.jl")
using .TRIPS
export
    calc_gravity,
    planet_structure,
    tidal_resp

include("thermal.jl")
using .ThE
export
    P1,
    one_layer_plnt,
    two_layer_plnt,
    init_profiles,
    temp,
    find_core

include("orbits.jl")
using .Orbit
export
    orbital_evolution,
    inclination_evolution

include("evo.jl")
using .Evolution
export
    planet_evo,
    tidal_evo

end # module
