#===============================================================================

    le_winding_generator.jl - LR, February 2021

    This is a script to create non-zero winding states from maximally flippable
    states in the [0,0,0] winding sector.

===============================================================================#
include("../src/typedefs.jl")
include("../src/io/io.jl")
include("../src/le_stuff/nonzero_wstate_generation.jl")

# Read params and make a logger.
param = read_config("wind_config.yml")
latt = LinkLattice(param["L"])

# Set the link type for the correct representation.
const LinkType = latt.S[end]*latt.d >= 63 ? LargeLinkState : SmallLinkState


# TODO: The base states need to be converted first. What happens for larger
# dataype - does it work natively? (the reading of the YAML file)
winding_states = ifncrease_winding(LinkType.(param["base_states"]), latt, ypos; increment=-1)

# Check output.
for state in winding_states
    println(Int(state))
end
