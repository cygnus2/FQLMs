#===============================================================================

    python_import.jl - LR, January 2021

    Some routines to read existant data produced with the Python version.

===============================================================================#
using HDF5
include("../typedefs.jl")


function read_lookup_tables(param::Dict{Any,Any})::Tuple{String,LookupDict,InvLookupDict}
    """ Reads the GL states for the specified winding sector(s).
    """
    latt = LinkLattice(param["L"])
    ws = haskey(param, "winding_sector") ? _winding_tag(param["winding_sector"], latt=latt) : "all-ws"

    # Read the correct states.
    states = param["low_energy_run"] ? _read_LE_states(param, ws) : _read_states(param, ws)

    # Construct the dictionaries.
    lookup_table = LookupDict()
    inverse_lookup_table = InvLookupDict()
    for k = 1:length(states)
        lookup_table[k] = states[k]
        inverse_lookup_table[states[k]] = k
    end

    return ws, lookup_table, inverse_lookup_table
end

function _read_states(param::Dict{Any,Any}, ws::String)::Array{LinkType,1}
    """ Reads full winding sectors - should only be called internally [from read_states()].
    """
    # Choose the correct file.
    filename = get(param, "state_file", nothing)
    if isnothing(filename)
        filename = param["working_directory"] * "/winding_states_" * _size_tag(param["L"]) * ".hdf5"
    end

    if ws == "all-ws"
        latt =  LinkLattice(param["L"])
        states = Array{LinkState,1}()
        for wsect in winding_sectors(latt)
            ws_states = h5open(filename, "r") do file
                read(file, _winding_tag(wsect))
            end
            states = vcat(states, ws_states)
        end
    else
        states = h5open(filename, "r") do file
            read(file, ws)
        end
    end
    return  Array{LinkType,1}(states)
end

# ---
# Reading and conversion of low-energy states.

# This was used in Python for the bitshfit - don't change.
const PY_bitshift = 63
function _convert_links_from_HDF5(data::Array{Integer,2})::Array{LargeLinkState}
    combined_data = LargeLinkState.(fill(0, size(data)[begin]))
    for (i, (x, y)) in enumerate(data)
        combined_data[i] = LargeLinkState(x)<<PY_bitshift + LargeLinkState(y)
    end
    return combined_data
end
_convert_links_from_HDF5(data::Array{Int64,1})::Array{SmallLinkState} = SmallLinkState.(data)


function _read_LE_states(param::Dict{Any,Any}, ws::String; combine=true)::Array{LinkType,1}
    """ Reads the GL states for the low-enegy sector up to a specified level. If
        the specified level is not found, the maximum will be used.

        Notes:
            - currently limited to the zero winding sector
            - should only be called internally [from read_states()]
    """
    # Choose the correct file.
    filename = get(param, "state_file", nothing)
    if isnothing(filename)
        filename = param["working_directory"] * "/le_states_" * _size_tag(param["L"]) * ".hdf5"
    end

    states = Dict()
    levels = Vector{Int}()
    if ws == "all-ws"
        h5open(filename, "r") do file
            for ds in file
                lv = parse(Int, HDF5.name(ds)[end])
                push!(levels, lv)
                states[lv] = sort(_convert_links_from_HDF5(ds[:]))
            end
        end
    else
        error("Only the zero-winding sector is implemented for LE!")
    end

    if maximum(levels) < param["maximum_excitation_level"]
        @warn "Using maximum excitation level that is lower than specified!" max_level = maximum(levels) specified=param["maximum_excitation_level"]
    end

    if combine
        all_states = Array{LinkType,1}()
        for (k, s) in states
            all_states = vcat(all_states, s)
        end
        return sort(all_states)
    else
        error("Separated LE states not yet implemented!")
    end
    return states
end
