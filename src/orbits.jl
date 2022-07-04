module Orbit

using Planets: Planet, Moon
using Structure: planet_mmotion
using Numerics: interpolate, rk4

export orbital_recession

function orbital_recession(plnt::Planet, mn::Moon, kl::AbstractVector,
                          Ql::AbstractVector, t::AbstractVector, t0::Real,
                          t1::Real; steps::Int=100)

    function dadt(t, a, p)

        k, dk = p.i(p.t, p.k, t)
        Q, dQ = p.i(p.t, p.Q, t)
        if(isnan(Q)) Q = 5e13 end

        a = a[1]
        n = p.n(p.GM, a)

        return 3*a * (abs(k)/Q) * n * (p.gm/p.GM) * (p.R/a)^5
    end

    param = (k=kl, Q=Ql, t=t, GM=plnt.GM, gm=mn.gm, R=plnt.R, i=interpolate,
             n=planet_mmotion)
    tspan = (t0, t1)
    sol = rk4(dadt, [mn.a], tspan, steps, p=param)

    return sol
end


end # module
