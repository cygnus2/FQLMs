#===============================================================================

    python_import.jl - LR, January 2021

    Some routines to read existant data produced with the Python version.

===============================================================================#
using HDF5
include("../typedefs.jl")

function _read_states(param, ws):: Array{LinkState,1}
    filename = param["working_directory"] * "/winding_states_" * _size_tag(param["L"]) * ".hdf5"
    if ws == "all-ws"
        latt =  LinkLattice(param["L"])
        states = Array{LinkState,1}()
        for wsect in winding_sectors(latt)
            ws_states = h5open(filename, "r") do file
                read(file, _winding_tag(collect(wsect)))
            end
            states = vcat(states, ws_states)
        end
    else
        states = h5open(filename, "r") do file
            read(file, ws)
        end
    end
    return  Array{LinkState,1}(states)
end

function read_lookup_tables(param)
    """ Reads the winding states
    """
    latt = LinkLattice(param["L"])
    ws = haskey(param, "winding_sector") ? _winding_tag(param["winding_sector"], latt=latt) : "all-ws"
    states = _read_states(param, ws)

    lookup_table = LookupDict()
    inverse_lookup_table = InvLookupDict()
    for k = 1:length(states)
        lookup_table[k] = states[k]
        inverse_lookup_table[states[k]] = k
    end

    return ws, lookup_table, inverse_lookup_table
end
