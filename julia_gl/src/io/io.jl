#===============================================================================

    io.jl - LR, January 2021

    Some I/O related stuff.

===============================================================================#
using YAML
include("../typedefs.jl")


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

function list_from_param(att::String, param::Dict{Any,Any})::Array{Number,1}
    """ Retrieves the list of parameters from the param dictionary.
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

function _size_tag(L)
    stag = ""
    for l in L
        stag *= "$(l)x"
    end
    return stag[begin:end-1]
end


function _winding_shift(latt::LinkLattice, ws::Array{Int,1})::Array{Int,1}
    """ Maps between the representation of winding numbers.
    """
    if latt.d == 2
        shift = reverse(latt.L) .÷ 2
    elseif latt.d == 3
        shift = [
            latt.L[2]*latt.L[3] ÷ 2;
            latt.L[1]*latt.L[3] ÷ 2;
            latt.L[1]*latt.L[2] ÷ 2
        ]
    else
        error("Dimension not implemented!")
    end
    return ws + shift
end


function _winding_tag(ws::Array{Int,1}; labels=['x', 'y', 'z'], latt::Union{LinkLattice,Nothing}=nothing)
    """ Returns the naming convention of the winding datasets.

        The shift is a lattice configuration L = [Lx,Ly,Lz], such
        that the HDF5 datasets may be resolved.
    """
    if !isnothing(latt)
        ws = _winding_shift(latt, ws)
    end
    wtag = ""
    for k=1:length(ws)
        wtag *= "w$(labels[k])_$(ws[k])-"
    end
    return wtag[begin:end-1]
end
