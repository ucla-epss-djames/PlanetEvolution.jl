module ThE
# [Th]ermal [E]volution

using Roots
using QuadGK
using DifferentialEquations

using Numerics
using Structure
using Planets
using Thermal

using ..PlanetEvolution: calc_gravity

export one_layer_plnt, two_layer_plnt, init_profiles, find_core

const P1 = 1 * bar_to_Pa

function one_layer_plnt(plnt::Planet, ρ::Function, T1::Real; t0::Real=0.0,
                        t1::Real=10.0)

    g, P = init_profiles(plnt, ρ)

    function dTdt(T, p, t)

        T_ef = p.T_ef(T, p.plnt.B)
        L_f = p.L(p.plnt.R, T_ef, p.plnt.T_eq)

        dT = L_f / p.I

        return dT
    end

    param = (plnt=plnt, ρ=ρ, P=P, i=interpolate)
    I = cons_mass(param, P1, 0, plnt.R) * -4*π * plnt.C_p

    u = T1
    tspan = (t0, t1*Gyr_to_sec)
    param = (plnt=plnt, I=I, T_ef=temp_effective, L=lumin_internal)

    prob = ODEProblem(dTdt, u, tspan, param)
    sol = solve(prob, reltol=1e-6, abstol=1e-7, Tsit5())

    T = sol.u
    t = sol.t
    t /= Gyr_to_sec

    return [t T]
end

function two_layer_plnt(plnt::Planet, ρ::Function, T1::Real; Ti::Real=0,
                        t0::Real=0.0, t1::Real=10.0)

    function dTdt(T, p, t)

        T1 = T[1]
        Ti = T[2]
        dT1 = dTi = 0

        T_ef = p.T_ef(T1, p.plnt.B)
        L = p.L(p.plnt.R, T_ef, p.plnt.T_eq)
        c, P_c, T_c, ρ_c, g_c = find_core(T1, p)

        if c == 0.0
            # if planet is still fluid

            I = p.I(p, P1, 0, p.plnt.R) * -4*π * p.plnt.C_p

            dT1 = L / I
            dTi = p.T(P_c, dT1, P1, p.plnt.∇)

        else
            # if planet has a core

            q = p.L_c(T1, Ti, c, P_c, T_c, ρ_c, g_c, P1, p.plnt)

            L_f = L - q.L_c

            I1 = p.I(p, P1, c, p.plnt.R) * -4*π * p.plnt.C_p
            I2 = p.I(p, P_c, 0, c) * -4*π * p.plnt.C_p

            dT1 = L_f / I1

            ∂c∂t = q.K * dT1

            dTi = q.L_c + 4*π*c^2*p.plnt.C_p*q.ΔT*∂c∂t*ρ_c + I2*Ti/T1*q.Γ*dT1
            dTi /= I2
        end

        return [dT1, dTi]
    end

    g, P = init_profiles(plnt, ρ)

    if(Ti == 0) Ti = temp_adiabat(P[1,2], T1, P1, plnt.∇) end

    u = [T1, Ti]
    tspan = (t0, t1*Gyr_to_sec)
    param = (plnt=plnt, ρ=ρ, g=g, P=P, find_core=find_core, I=cons_mass,
             T_ef=temp_effective, T=temp_adiabat, T_m=temp_melting,
             L=lumin_internal, L_c=lumin_core, i=interpolate)

    prob = ODEProblem(dTdt, u, tspan, param)
    sol = solve(prob, reltol=1e-7, abstol=1e-8, Tsit5())
    T1 = sol[1,:]
    Ti = sol[2,:]
    t = sol.t
    t /= Gyr_to_sec

    return [t T1 Ti]
end

function init_profiles(plnt::Planet, ρ::Function)

    # radial profile
    r = 0.01:100:plnt.R
    l = length(r)

    # gravity profile
    g = zeros(l)
    mass = 0
    r0 = 0
    for i in 1:l

        r1 = r[i]
        if(i != 1) r0 = r[i-1] end
        mass, gi = calc_gravity(r0, r1, mass, ρ)
        g[i] = gi

    end

    # pressure profile
    dP(u, p, x) = dPdr(ρ(x), interpolate(r, g, x)[1])
    rspan = (plnt.R, 0)
    u0 = 0

    prob = ODEProblem(dP, u0, rspan)
    sol = solve(prob, reltol=1e-8, abstol=1e-10, Tsit5())
    P = sol.u[end:-1:1]
    x = sol.t[end:-1:1]

    return ([r g], [x P])
end

function cons_mass(p, P::Real, a::Real, b::Real)

    A, err = quadgk(r -> layer_density(r, p.i(p.P[:,1], p.P[:,2], r)[1], P,
                                       p.plnt.∇, p.ρ(r)), a, b)

    return A
end

function find_core(T1::Real, p)

    # cross section of melting and adiabat
    f(P) = p.T(P, T1, P1, p.plnt.∇) - p.T_m(P, p.plnt.P0, p.plnt.T0, p.plnt.a,
                                            p.plnt.b)

    # finding cross intersection of melting and adiabat
    x = find_zeros(f, p.plnt.P0, p.P[1,2])

    if isempty(x)
        # if no cross section set core pressure to central pressure
        P_c = p.P[1,2]
        c = 0.0
    else
        P_c = x[1]

        # use core pressure to interpolate the radius of the core
        c, dc = interpolate(p.P[:,2], p.P[:,1], P_c)

    end

    # gather the rest of the surface core values
    T_c = p.T(P_c, T1, P1, p.plnt.∇)
    ρ_c = p.ρ(c)
    g_c = interpolate(p.g[:,1], p.g[:,2], c)[1]

    return (c, P_c, T_c, ρ_c, g_c)
end

end # module
