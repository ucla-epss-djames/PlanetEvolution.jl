using Tidal
using Planets

export planet_structure, tidal_resp

"""
    planet_structure(plnt::Planet, mn::Moon, data::Matrix)

Generates structure matrix of a planet for tidal calculation. The matrix is a 4
column matrix containing radius, complex shear modulus, gravity, and density.
"""
function planet_structure(plnt::Planet, mn::Moon, data::Matrix)

    mass = 0.0
    r0 = 0.0
    layers = plnt.layers
    ω = plnt.ω
    μ_f = plnt.μ_f
    model = plnt.rhea_model

    n = planet_mmotion(plnt.GM, mn.a)
    χ = 2 * abs(ω - n)

    sd = zeros(Complex, layers, 4)

    for i in 1:layers

        r1 = data[i,1]
        ρ = data[i,2]
        μ = data[i,3]
        η = data[i,4]

        if(i != 1) r0 = real(sd[i-1,1]) end
        mass, g = calc_gravity(r0, r1, mass, ρ)

        sd[i,4] = ρ
        sd[i,3] = g
        sd[i,2] = planet_cmu(μ, χ, η, r1, g, ρ, μ_f, model)
        sd[i,1] = r1

    end

    return sd, mass
end

"""
    tidal_resp(plnt::Planet, mn::Moon, data::Matrix, flag::Bool; l::Int=2)

Calculates the tidal response of a planet.
"""
function tidal_resp(plnt::Planet, mn::Moon, data::Matrix, flag::Bool; l::Int=2)

    sd, mass = planet_structure(plnt, mn, data)

    normalize!(mass, real(sd[end,1]), sd)
    tidal = propagator_method(l, plnt.layers, sd, flag)

end
