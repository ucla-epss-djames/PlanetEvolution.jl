using Roots
using QuadGK
using OrdinaryDiffEq
using Interpolations

using Numerics
using Planets: Planet, calc_gravity, thermal_inertia, temp_effective, lumin_internal


export P1
export one_layer_plnt, two_layer_plnt, init_profiles, temp, find_core

const P1 = 1 * bar_to_Pa

"""
    one_layer_plnt(plnt::Planet, ρ::Function, T1::Real; t0::Real=0.0,
                   t1::Real=10.0)

Solves out thermal evolution for a fluid planet. Returns a two column matrix where
the first column is time [Gyr] and the second column is temperature [K].

# Arguments
- `plnt::Planet` - planet parameters
- `ρ::Function`  - density profile
- `T1::Real`     - initial temperature for solver
- `t0::Real`     - time start of thermal evolution
- `t1::Real`     - time end of thermal evolution
"""
function one_layer_plnt(plnt::Planet, ρ::Function, T1::Real; t0::Real=0.0,
                        t1::Real=10.0, reltol::Real=1e-6, abstol::Real=1e-7)

    g, P = init_profiles(plnt, ρ)

    function dTdt(T, p, t)

        T_ef = p.T_ef(T, p.plnt.B)
        L_f = p.L(p.plnt.R, T_ef, p.plnt.T_eq)

        dT = L_f / p.I

        return dT
    end

    param = (plnt=plnt, ρ=ρ, P=P, i=Numerics.interpolate)
    I = thermal_inertia(param, P1, 0, plnt.R) * -4*π * plnt.C_p

    u = T1
    tspan = (t0, t1*Gyr_to_sec)
    param = (plnt=plnt, I=I, T_ef=temp_effective, L=lumin_internal)

    prob = ODEProblem(dTdt, u, tspan, param)
    sol = solve(prob, reltol=reltol, abstol=abstol, Tsit5())

    T = sol.u
    t = sol.t
    t /= Gyr_to_sec

    return [t T]
end

"""
    two_layer_plnt(plnt::Planet, ρ::Function, T1::Real; Ti::Real=0, t0::Real=0.0,
                   t1::Real=10.0)

Solves out thermal evolution for a two layer planet. Returns a three column matrix
where the first column is time [Gyr], the second is temperature of the planet [K],
and the thrid is temperature of the thermal boundary layer [K].

# Arguments
- `plnt::Planet` - planet parameters
- `ρ::Function`  - density profile
- `T1::Real`     - initial temperature for solver
- `Ti::Real`     - initial temperature for thermal booundary layer
- `t0::Real`     - time start of thermal evolution
- `t1::Real`     - time end of thermal evolution
- `reltol::Real` - relative tolerance
- `abstol::Real` - absolute tolerance
"""
function two_layer_plnt(plnt::Planet, ρ::Function, T1::Real; Ti::Real=0,
                        t0::Real=0.0, t1::Real=10.0, reltol::Real=1e-7,
                        abstol::Real=1e-8)

    function dTdt!(dT, T, p, t)

        T_ef = p.T_ef(T[1], p.plnt.B)
        L = p.L(p.plnt.R, T_ef, p.plnt.T_eq)
        c, P_c, T_c, ρ_c, g_c = find_core(T[1], p)

        if c == 0.0
            # if planet is still fluid

            I = p.I(p, P1, 0, p.plnt.R) * -4*π * p.plnt.C_p

            dT[1] = L / I
            dT[2] = p.T(P_c, dT[1], P1, p.plnt.∇)

        else
            # if planet has a core

            q = p.L_c(T[1], T[2], c, P_c, T_c, ρ_c, g_c, P1, p.plnt)

            L_f = L - q.L_c

            I1 = p.I(p, P1, c, p.plnt.R) * -4*π * p.plnt.C_p
            I2 = p.I(p, P_c, 0, c) * -4*π * p.plnt.C_p

            dT[1] = L_f / I1

            ∂c∂t = q.K * dT[1]

            dT[2] = q.L_c + 4*π*c^2*p.plnt.C_p*q.ΔT*∂c∂t*ρ_c +
                I2*T[2]/T[1]*q.Γ*dT[1]
            dT[2] /= I2
        end

        nothing
    end

    g, P = init_profiles(plnt, ρ)

    if(Ti == 0) Ti = temp_adiabat(P[1,2], T1, P1, plnt.∇) end

    u = [T1, Ti]
    tspan = (t0, t1*Gyr_to_sec)
    param = (plnt=plnt, ρ=ρ, g=g, P=P, find_core=find_core, I=thermal_inertia,
             T_ef=temp_effective, T=temp_adiabat, T_m=temp_melting,
             L=lumin_internal, L_c=lumin_core, i=Numerics.interpolate)

    prob = ODEProblem(dTdt!, u, tspan, param)
    sol = solve(prob, reltol=reltol, abstol=abstol, Vern6())
    T1 = sol[1,:]
    Ti = sol[2,:]
    t = sol.t
    t /= Gyr_to_sec

    return [t T1 Ti]
end

"""
    init_profiles(plnt::Planet, ρ::Function)

Generates a gravity and pressure profile for a planet.

# Arguments
- `plnt::Planet` - planet parameters
- `ρ::Function`  - density profile
"""
function init_profiles(plnt::Planet, ρ::Function)
    # updating mass and returning a new plnt

    # radial profile
    r = 0.01:plnt.layers:plnt.R
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
    g = CubicSplineInterpolation(r, g, extrapolation_bc=Line())
    println(plnt.name, " mass: ", mass, " kg")

    # pressure profile
    dP(u, p, x) = dPdr(ρ(x), g(x))
    rspan = (plnt.R, 0)
    u0 = 0

    prob = ODEProblem(dP, u0, rspan)
    sol = solve(prob, reltol=1e-8, abstol=1e-10, OwrenZen3())
    P = sol.u[end:-1:1]
    x = sol.t[end:-1:1]

    return (g, [x P])
end

"""
    temp(r::Real, T1::Real, Ti::Real, P_c::Real, p)

Calculates the adiabatic temperature based on if `r` is inside the core or the
envelope. Refer to Stixrude et al. 2021 (eq 12).

# Arguments
- `r::Real`   - radius
- `T1::Real`  - temperature
- `Ti::Real`  - thermal boundary temperature
- `P_c::Real` - pressure at top of the core
- `p`         - tuple holding several parameters
"""
function temp(r::Real, T1::Real, Ti::Real, P_c::Real, p)
    P, dP = Numerics.interpolate(p.P[:,1], p.P[:,2], r)
    if P < P_c
        return temp_adiabat(P, T1, P1, p.plnt.∇)
    else
        return temp_adiabat(P, Ti, P_c, p.plnt.∇)
    end
end

"""
    thermal_inertia(p, P::Real, a::Real, b::Real)

Integrates the `layer_density` for the luminosity of a planet.

# Arguments
- `p`        - tuple holding several parameters
- `P::Real`  - pressure
- `r0::Real` - radial start
- `r1::Real` - radial end
"""
function thermal_inertia(p, P::Real, r0::Real, r1::Real)

    A, err = quadgk(r -> layer_density(r, p.i(p.P[:,1], p.P[:,2], r)[1], P,
                                       p.plnt.∇, p.ρ(r)), r0, r1)

    return A
end

"""
    find_core(T1::Real, p)

Finds top of core values based on intersection of adiabat and melting temperature
curves.

# Arguments
- `T1::Real` - temperature
- `p`        - tuple holding several parameters
"""
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
        c, dc = p.i(p.P[:,2], p.P[:,1], P_c)

    end

    # gather the rest of the surface core values
    T_c = p.T(P_c, T1, P1, p.plnt.∇)
    ρ_c = p.ρ(c)
    g_c = p.g(c)

    return (c, P_c, T_c, ρ_c, g_c)
end
