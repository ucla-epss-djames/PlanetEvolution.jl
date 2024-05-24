using StructTypes
using JSON3
using Planets

export write_json, read_json

struct Data_Schema

    # meta
    date    ::    String

    # planetary parameters
    plnt       ::    Planet
    mn         ::    Moon
    density    ::    String

    # evolutionary output data
    cols_one_lyr    ::    Vector{String}
    size_one_lyr    ::    Vector{Int64}
    data_one_lyr    ::    Matrix{Float64}
    cols_two_lyr    ::    Vector{String}
    size_two_lyr    ::    Vector{Int64}
    data_two_lyr    ::    Matrix{Float64}

end
StructTypes.StructType(::Type{Data_Schema}) = StructTypes.Struct()

function write_json(d, path)

    open(path * "planet_evo-" * d["date"] * ".json", "w") do f
        JSON3.write(f, d, allow_inf=true)
        println(f)
    end

end

function read_json(path)

    d = JSON3.read(path, allow_inf=true)

    data_one_lyr = convert(Matrix{Float64},
                           reshape(d.data_one_lyr, (d.size_one_lyr...)))
    data_two_lyr = convert(Matrix{Float64},
                           reshape(d.data_two_lyr, (d.size_two_lyr...)))

    plnt = d.plnt
    mn = d.mn
    return Data_Schema(d.date, Planet(plnt.name, plnt.layers, plnt.R, plnt.M,
                                      plnt.GM, plnt.ω, plnt.C_p, plnt.α, plnt.k,
                                      plnt.κ, plnt.T0, plnt.P0, plnt.a, plnt.b,
                                      plnt.B, plnt.∇, plnt.T_eq, plnt.T_ef, plnt.T1,
                                      plnt.η0, plnt.A, plnt.Ra,
                                      (plnt.μ_f[1], plnt.μ_f[2], plnt.μ_f[3]),
                                      plnt.rhea_model),
                       Moon(mn.name, mn.m, mn.gm, mn.a, mn.i, mn.retro),
                       d.density, d.cols_one_lyr, d.size_one_lyr, data_one_lyr,
                       d.cols_two_lyr, d.size_two_lyr, data_two_lyr)

end
