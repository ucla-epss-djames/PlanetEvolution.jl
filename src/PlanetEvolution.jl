module PlanetEvolution

greet() = print("Hello World!")

include("tidal.jl")
using .TRIPS

include("thermal.jl")
using .OTTER


end # module
