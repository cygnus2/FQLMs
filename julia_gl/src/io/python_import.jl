#===============================================================================

    python_import.jl - LR, January 2021

    Some routines to read existant data produced with the Python version.

===============================================================================#
using HDF5
include("../typedefs.jl")


function read_lookup_tables(param::Dict{Any,Any})::Tuple{LookupDict,InvLookupDict}
    """ Reads the GL states for the specified winding sector(s).
    """
    # Read the correct states.
    states = param["low_energy_run"] ? _read_LE_states(param) : _read_states(param)

    # Construct the dictionaries.
    lookup_table = LookupDict()
    inverse_lookup_table = InvLookupDict()
    for k = 1:length(states)
        lookup_table[k] = states[k]
        inverse_lookup_table[states[k]] = k
    end

    return lookup_table, inverse_lookup_table
end

function _read_states(param::Dict{Any,Any})::Array{LinkType,1}
    """ Reads full winding sectors - should only be called internally [from read_states()].
    """
    if param["ws_label"] == "all-ws"
        latt =  LinkLattice(param["L"])
        states = Array{LinkState,1}()
        for wsect in winding_sectors(latt)
            ws_states = h5open(param["state_file"], "r") do file
                read(file, _winding_tag(wsect))
            end
            states = vcat(states, ws_states)
        end
    else
        states = h5open(param["state_file"], "r") do file
            read(file, param["ws_label"])
        end
    end
    return  Array{LinkType,1}(states)
end

# ---
# Reading and conversion of low-energy states.

# This was used in Python for the bitshfit - don't change.
const PY_bitshift = 63
function _convert_links_from_HDF5(data::Array{Int64,2})::Array{LargeLinkState}
    combined_data = LargeLinkState.(fill(0, size(data)[end]))
    for k = 1:length(combined_data)
        (x, y) = data[:,k]
        combined_data[k] = LargeLinkState(x)<<PY_bitshift + LargeLinkState(y)
    end
    return combined_data
end
_convert_links_from_HDF5(data::Array{Int64,1})::Array{SmallLinkState} = SmallLinkState.(data)


function _read_LE_states(param::Dict{Any,Any}; combine=true)::Array{LinkType,1}
    """ Reads the GL states for the low-enegy sector up to a specified level. If
        the specified level is not found, the maximum will be used.

        Notes:
            - currently limited to the zero winding sector
            - should only be called internally [from read_states()]
    """
    # Choose the correct file.
    states = Dict()
    levels = Vector{Int}()
    if param["ws_label"] == "all-ws"
        h5open(param["state_file"], "r") do file
            for ds in file
                lv = parse(Int, HDF5.name(ds)[end])
                if lv <= param["maximum_excitation_level"]
                    push!(levels, lv)
                    # Julia lacks the ellipsis feature, therfore here goes a slight
                    # misuse of mutiple dispatch.
                    if length(ds)[end] > 0
                        if length(size(ds)) == 2
                            states[lv] = sort(_convert_links_from_HDF5(ds[:,:]))
                        else
                            states[lv] = sort(_convert_links_from_HDF5(ds[:]))
                        end
                    end
                end
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
