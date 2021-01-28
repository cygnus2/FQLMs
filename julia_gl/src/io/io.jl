#===============================================================================

    io.jl - LR, January 2021

    Some I/O related stuff.

===============================================================================#
using YAML
using Printf
include("../typedefs.jl")

function read_config(filename::String)
    conf = YAML.load(open(filename))
    return conf
end

function list_from_param(att::String, param::Dict{Any,Any})::Array{Number,1}
    """ Retrieves the list of parameters from the param dictionary (or throws
        an error, if not possible).
    """
    if haskey(param, att*"_range")
        (b, e, len) = param[att*"_range"]
        return range(b, e; length=len)
    elseif haskey(param, att*"_list")
        return param[att*"_list"]
    else
        return [param[att]]
    end
end

function _lambda_tag(lambda::Number)::String
    """ Central place that toggls how the datasets are named in the HDF5 output.
    """
    return "_lam_"*(@sprintf "%.6f" lambda)
end
