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


function _convert

function _read_states(param::Dict{Any,Any}, ws::String)::Array{LinkState,1}
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
    return  Array{LinkState,1}(states)
end


# function _combine(data; bit_shift=63):
#     """ Converter method required to map back two integers with maximum
#         64 bit (like HDF5 has it) to a larger datatype. This is needed to store
#         the states for lattices with more than 63 spins.
#     """
#     combined_data =
#     combined_data = [2**70]*len(data)
#     for i, (x, y) in enumerate(data):
#         combined_data[i] = (int(x)<<bit_shift) + int(y)
#     return sorted(combined_data)
# end
#
#
#
# function _read_LE_states(param::Dict{Any,Any}, ws::String; combine=true):: Array{LinkState,1}
#     """ Reads the GL states for the low-enegy sector up to a specified level. If
#         the specified level is not found, the maximum will be used.
#
#         Notes:
#             - currently limited to the zero winding sector
#             - should only be called internally [from read_states()]
#     """
#     # Choose the correct file.
#     filename = get(param, "state_file", nothing)
#     if isnothing(filename)
#         filename = param["working_directory"] * "/le_states_" * _size_tag(param["L"]) * ".hdf5"
#     end
#
#         states = Dict()
#         if ws == "all-ws"
#             h5open(filename, "r") do file
#                 for ds in file
#                     lv = parse(Int, HDF5.name(ds)[end])
#                     states[lv] = ds[:]
#                 end
#             end
#         else
#             error("Only the zero-winding sector is implemented for LE!")
#         end
#
#         states = {}
#         with hdf.File(state_file, 'r') as f:
#             levels = []
#             for g in f:
#                 l = int(g.split("_")[-1])
#                 levels.append(l)
#                 if l <= max_level:
#                     s = f[g][...]
#                     if len(s.shape) == 2:
#                         states[l] = LowEnergyGLSimulation._combine(s)
#                     else:
#                         states[l] = sorted(map(int, s))
#
#             if sorted(levels)[-1] < max_level:
#                 self.log(f"Warning: maximal level not reachable from list! Taking maximum of {l}.")
#                 new_max = l
#             else:
#                 new_max = max_level
#
#
#         if combine
#             return np.concatenate([states[k] for k in sorted(states.keys()) if len(states[k])])
#         else
#             error("Separated LE states not yet implemented!")
#         end
#         return states
#
#     except OSError:
#         self.log(f'Could not find file {state_file}')
#         return [], 0
#
# end
#
#
# def _get_hamiltonian_file(self, prefactor='LE_hamiltonian', default=None):
#     return GLSimulation._get_hamiltonian_file(self, prefactor=prefactor, default=default)
#
# def _get_result_file(self, prefactor='LE_results', default=None):
#     return GLSimulation._get_result_file(self, prefactor=prefactor, default=default)
#
# def _get_state_file(self, prefactor='le_states', default=None):
#     return GLSimulation._get_state_file(self, prefactor=prefactor, default=default)
