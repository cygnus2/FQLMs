#===============================================================================

    opeator_io.jl - LR, January 2021

    Read operator in a safe way.

===============================================================================#
using YAML
include("../typedefs.jl")


function read_operators(filename="operator_masks.yml")::Dict{String,Dict{String,Operator}}
    """ Retrieves the pre-stored operators from the Python generated YAML file.
    """
    raw_operators = YAML.load(open(filename))

    op_dict = Dict{String,Dict{String,Operator}}()
    for (latt, raw_latt_ops) in raw_operators
        latt_ops = Dict{String,Operator}()
        for (name, bit_map) in raw_latt_ops
            latt_ops[name] = Operator(name, bit_map)
        end
        op_dict[latt] = latt_ops
    end
    return op_dict
end
