#===============================================================================

    python_import.jl - LR, January 2021

    Some routines to read existant data produced with the Python version.

===============================================================================#
using HDF5
include("io.jl")


function read_states(param, ws)
    filename = param["working_directory"] * "/winding_states_" * _size_tag(param["L"]) * ".hdf5"

    if ws == "all-ws"
        return nothing
    else
        states = h5open(filename, "r") do file
            read(file, ws)
        end
    end
    return states
end

    #
    #
    # def read_states(self, file=None, merged=True):
    #     """ Read the GLS from file.
    #     """
    #     state_file = self._get_state_file(default=file)
    #     try:
    #         states = []
    #         if self.ws == 'all-ws':
    #             with hdf.File(state_file, 'r') as f:
    #                 for ws in f:
    #                     if merged:
    #                         states += list(f[ws][...])
    #                     else:
    #                         states.append([ws, list(f[ws][...])])
    #         else:
    #             with hdf.File(state_file, 'r') as f:
    #                 states = f[self.ws][...]
    #
    #         self.log(f'Read Fock states from {state_file}')
    #         return states
    #
    #     except (OSError):
    #         self.log(f'Could not find file {state_file}')
    #         return []
