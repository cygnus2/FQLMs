#===============================================================================

    plaquette.jl - LR, Jan 2021

    Representation of plaquettes with an immutable struct.

    The convention is the following: link link_1 and link_2 as well as link_3
    and link_4 are connected through the hopping like

                                  link3
                                O--------O
                                |        |
                          link4 |        | link2
                                |        |
                                O--------O
                                  link1

        (ASCII art shamelessly stolen from Debasish)

    This list will later be used when we loop through the states to
    construct the Hamiltonian.

===============================================================================#
struct Plaquette
    links::Vector{LinkIndex} # Holds the four participating link indicies.
    mask::LinkState # Represents the links as a LinkState.

    Plaquette(links::Vector{LinkIndex}) = begin
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
Base.length(p::Plaquette) = 4
function Base.iterate(p::Plaquette, state::Integer=1)
    state > length(p) ? nothing : (p.links[state], state+1)
end

# Some convenient stuff.
Base.maximum(p::Plaquette) = maximum(p.links)
