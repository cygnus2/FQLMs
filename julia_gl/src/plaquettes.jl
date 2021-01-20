#===============================================================================

    rep_struct.jl - LR, Jan 2021

    Representation of plaquettes with an immutable struct.

===============================================================================#
include("../state_reps/rep_integer.jl")

struct Plaquette
    links::Vector{UInt8} # Holds the four participating link indicies.
    mask::LinkState # Represents the links as a LinkState.

    Plaquette(links::Vector{UInt8}) = begin
        mask = LinkState(0)
        for l in links
            mask = create(mask, l)
            if isnothing(mask)
                error("Doubly provided link in plaquette creation.")
            end
        end
        new(links, mask)
    end
end

# For indexing.
function Base.getindex(p::Plaquette, i::Integer)
    return p.links[i]
end
Base.firstindex(p::Plaquette) =  1
Base.lastindex(p::Plaquette) = 4

# For iteration.
Base.length(S::Plaquette) = 4
function Base.iterate(p::Plaquette, state::Integer=1)
    state > length(p) ? nothing : (p.links[state], state+1)
end
