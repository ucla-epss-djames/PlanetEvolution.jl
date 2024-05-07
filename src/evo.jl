module Evolution

using Printf

using Numerics
using Thermal
using Planets

using ..PlanetEvolution: one_layer_plnt, two_layer_plnt, init_profiles, find_core, P1, tidal_resp, write_json

export planet_evo, tidal_evo

function planet_evo(plnt, rho, rho_path, T1, mn, path)

    d = Dict("date" => get_date(), "plnt" => plnt, "mn" => mn,
             "density" => rho_path,
             "cols_one_lyr" => ["time [s]", "effective temp [K]"])

    println("RUNNING ONE LAYER")
    one_lyr = one_layer_plnt(plnt, rho, T1)
    l = length(one_lyr[:,1])
    d["size_one_lyr"] = [l,2]
    d["data_one_lyr"] = one_lyr

    println("RUNNING TWO LAYER MODEL")
    two_lyr = two_layer_plnt(plnt, rho, T1)
    l = length(two_lyr[:,2])
    d["size_two_lyr"] = [l,13]
    d["cols_two_lyr"] = union(d["cols_one_lyr"], ["internal temp [K]",
                                                  "surface core temp [K]",
                                                  "TBL temp [K]",
                                                  "surface core pres [Pa]",
                                                  "core radius [m]",
                                                  "surface core rho [kg m⁻³]",
                                                  "surface core grav",
                                                  "Rayleigh Number [-]",
                                                  "TBL size [m]",
                                                  "Tidal Love Number k₂ [-]",
                                                  "Tidal Quality [-]"])

    println("INITIALIZING PROFILES")
    g, P = init_profiles(plnt, rho)
    p = (plnt=plnt, ρ=rho, g=g, P=P, T=temp_adiabat, T_m=temp_melting,
         i=Numerics.interpolate)

    vals = [find_core(T, p) for T in two_lyr[:,2]]
    c = [val[1] for val in vals]
    P_c = [val[2] for val in vals]
    T_c = [val[3] for val in vals]
    DT = two_lyr[:,3] .- T_c
    rho_c = [val[4] for val in vals]
    g_c = [val[5] for val in vals]

    kl, Ql = tidal_evo(plnt, rho, two_lyr, mn, p, c)

    l = length(two_lyr[:,1])
    s = findnext(x -> x > 0, c, 1)
    q = [lumin_core(two_lyr[i,2], two_lyr[i,3], c[i], P_c[i], T_c[i], rho_c[i],
                    g_c[i], P1, plnt) for i in s:l]
    Ra = [(i < s ? 0. : q[i-s+1].Ra) for i in 1:l]
    δ = c .* (plnt.Ra ./ Ra) .^ (1/3)

    d["data_two_lyr"] = [two_lyr T_c DT P_c c rho_c g_c Ra δ kl Ql]

    write_json(d, path)

end

function tidal_evo(plnt, rho, T2, mn, p, c)

    mu(T,P) = (101 + P*1e-9/2.41 - (T - 1650)/23.4) * 1e9

    l = findnext(x -> x >= 10, T2[:,1], 1)

    println(plnt.name, " EVOLVING WITH ", mn.name)

    ri = 0.1
    println("PROCESSING TIDAL EVOLUTION")
    output = zeros(l,2)
    for i in 1:l

        plntt = zeros(p.plnt.layers, 4)
        layers_c = floor(Int, (c[i] / p.plnt.R) * p.plnt.layers)
        if(layers_c == 0 || layers_c <= 100)
            plntt[:,1] = range(ri, p.plnt.R, length=p.plnt.layers)
        else
            layers_m = p.plnt.layers - layers_c
            plntt[1:layers_c,1] = range(ri, c[i], length=layers_c)
            dr = plntt[2,1] - plntt[1,1]
            plntt[layers_c+1:p.plnt.layers,1] = range(c[i] + dr, p.plnt.R,
                                                     length=layers_m)
            η = zeros(layers_c)
            μ = zeros(layers_c)
            for j in 1:layers_c
                P_r, dP = p.i(p.P[:,1], p.P[:,2], plntt[j,1])
                T_m = temp_melting(P_r, p.plnt.P0, p.plnt.T0, p.plnt.a, p.plnt.b)
                T = temp_adiabat(P_r, T2[i,2], P1, p.plnt.∇)
                η[j] = planet_eta(p.plnt.η0, p.plnt.A, T_m, T)
                μ[j] = mu(T, P_r)
            end

            plntt[1:layers_c,4] = η
            plntt[1:layers_c,3] = μ
        end

        plntt[:,2] = p.ρ.(plntt[:,1])

        tidal = tidal_resp(p.plnt, mn, plntt, false)
        output[i,1] = real(tidal[end,1])
        output[i,2] = imag(tidal[end,1])
    end


    kl = output[1:l,1]
    Ql = abs.(output[1:l,1] ./ output[1:l,2])

    return kl, Ql
end


end # module
