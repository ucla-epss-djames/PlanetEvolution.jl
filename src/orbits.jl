module Orbit

using Planets: Planet, Moon
using Structure: planet_mmotion
using Numerics: interpolate, rk4

export orbital_evolution, inclination_evolution

"""
    orbital_evolution(plnt::Planet, mn::Moon, kl::AbstractVector,
                      Ql::AbstractVector, t::AbstractVector, t0::Real,
                      t1::Real; steps::Int=100)

Calculates the orbital evolution of a planet's moon and describes the tides
raised on the planet. Refer to Murray's Solar System Dynamics Eq. 4.160.

# Arguments
- `plnt::Planet`       - planet parameters
- `mn::Moon`           - moon parameters
- `kl::AbstractVector` - kl tidal number evolution
- `Ql::AbstractVector` - Ql tidal quality evolution
- `t::AbstractVector`  - time span of `kl` and `Ql`
- `t0::Real`           - time start of evolution
- `t1::Real`           - time end of evolution
- `steps::Int=100`     - number of steps to take in solver
"""
function orbital_evolution(plnt::Planet, mn::Moon, kl::AbstractVector,
                           Ql::AbstractVector, t::AbstractVector, t0::Real,
                           t1::Real; steps::Int=100)

    function dadt(t, a, p)

        k, dk = p.i(p.t, p.k, t)
        Q, dQ = p.i(p.t, p.Q, t)
        if(isnan(Q)) Q = 5e13 end

        a = a[1]
        n = p.n(p.GM, a)

        if(p.retro)
            # if moon has orbital retrograde
            return -3*a * (abs(k)/Q) * n * (p.gm/p.GM) * (p.R/a)^5
        else
            return sign(p.ω - n) *3*a * (abs(k)/Q) * n * (p.gm/p.GM) * (p.R/a)^5
        end

    end

    param = (k=kl, Q=Ql, t=t, GM=plnt.GM, gm=mn.gm, R=plnt.R, ω=plnt.ω,
             retro=mn.retro, i=interpolate, n=planet_mmotion)
    tspan = (t0, t1)
    sol = rk4(dadt, [mn.a], tspan, steps, p=param)

    return sol
end

"""
    inclination_evolution(plnt::Planet, mn::Moon, kl::AbstractVector,
                          Ql::AbstractVector, t::AbstractVector, t0::Real,
                          t1::Real; steps::Int=100)

Calculates the inclination evolution of a planet's moon. Refer to
Chyba et al. '89, eq. 11, for more description.

# Arguments
- `plnt::Planet`       - planet parameters
- `mn::Moon`           - moon parameters
- `kl::AbstractVector` - kl tidal number evolution
- `Ql::AbstractVector` - Ql tidal quality evolution
- `t::AbstractVector`  - time span of `kl` and `Ql`
- `t0::Real`           - time start of evolution
- `t1::Real`           - time end of evolution
- `steps::Int=100`     - number of steps to take in solver
"""
function inclination_evolution(plnt::Planet, mn::Moon, kl::AbstractVector,
                               Ql::AbstractVector, t::AbstractVector, t0::Real,
                               t1::Real; steps::Int=100)

    function didt(t, a, p)

        k, dk = p.i(p.t, p.k, t)
        Q, dQ = p.i(p.t, p.Q, t)
        if(isnan(Q)) Q = 5e13 end

        a = a[1]
        n = p.n(p.GM, a)
        sin_i = (p.a/a)^(1/4) * sin(p.i0)

        return -3/4 * sin_i * (abs(k)/Q) * n * (p.gm/p.GM) * (p.R/a)^5
    end

    param = (k=kl, Q=Ql, t=t, GM=plnt.GM, gm=mn.gm, R=plnt.R, a=mn.a, i0=mn.i,
             i=interpolate, n=planet_mmotion)
    tspan = (t0, t1)
    sol = rk4(didt, [mn.a], tspan, steps, p=param)

    return sol
end


end # module
