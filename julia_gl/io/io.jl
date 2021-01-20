#===============================================================================

    io.jl - LR, January 2021

    Some I/O related stuff.

===============================================================================#
using YAML

function read_config(filename::String)
    conf = YAML.load(open(filename))
    if _sanity_checks!(conf)
        return conf
    end
    return nothing
end

function _sanity_checks!(conf)::Bool
    if !haskey(conf, "working_directory")
        conf["working_directory"] = "./"
    end
    return true
end

function _size_tag(L)
    stag = ""
    for l in L
        stag *= "$(l)x"
    end
    return stag[begin:end-1]
end
