module Util

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
    size_one_lyr    ::    Vector{Float64}
    data_one_lyr    ::    Vector{Vector{Float64}}
    cols_two_lyr    ::    Vector{String}
    size_two_lyr    ::    Vector{Float64}
    data_two_lyr    ::    Vector{Vector{Float64}}

end
StructTypes.StructType(::Type{Data_Schema}) = StructTypes.Struct()

function write_json(d, path)

    open(path * "planet_evo-" * d["date"] * ".json", "w") do f
        JSON3.write(f, d, allow_inf=true)
        println(f)
    end

end

function read_json(path)

    d = JSON3.read(path, Data_Schema, allow_inf=true)

    d.data_one_lyr = convert(Matrix{Float64},
                             reshape(d.data_one_lyr, (d.size_one_lyr)))
    d.data_two_lyr = convert(Matrix{Float64},
                             reshape(d.data_two_lyr, (d.size_two_lyr)))

    return d
end

end # module
