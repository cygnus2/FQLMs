#===============================================================================

    _set_type.jl - LR, May 2021

    Fixes the datatype for the Lattice representation.

===============================================================================#
include("../../src/typedefs.jl")


function set_type(param::Dict{Any,Any})
    # Create the lattice structure that we're working on.
    local latt = LinkLattice(param["L"])

    # (when dealing with links, this type should be used to convert the datatype to
    # be consistent - otherwise there'll likely be an error or *very* subtle bugs)
    LinkType = latt.S[end]*latt.d >= 63 ? LargeLinkState : SmallLinkState
    @info "Set the datatype for link representation" type=LinkType
    return LinkType
end
