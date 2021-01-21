#===============================================================================

    state.jl - LR, Jan 2021

    Representation of states is done via an Integer.

===============================================================================#
# That's the occupation state of a lattice of links. Currently limited to
# 64 links, but we can have more with larger datatypes.
const LinkState = UInt64


function create(s::LinkState, i::Integer)::Union{Nothing,LinkState}
    """ Destroys a particle at position i, false if there's a particle already.
    """
    shift = 1 << (i-1)
    if s & shift == 0
        return s ⊻ shift
    else
        return nothing
    end
end

function annihilate(s::LinkState, i::Integer)::Union{Nothing,LinkState}
    """ Destroys a particle at position i, false if there's no particle to destroy.
    """
    shift = 1 << (i-1)
    if s & shift > 0
        return s ⊻ shift
    else
        return nothing
    end
end

function count_occupancies(s::LinkState, src::LinkIndex, dst::LinkIndex)
    """ Sums all occupancies between index src and dst, both inclusive.
        Assumes dst >= src.
    """
    # return sum([(s >>> k) & 1 for k=0:63]) --> slower!
    c = 0
    for k=src-1:dst-1
        c += (s >>> k) & 1
    end
    return c
end


function sign(s::LinkState, src::Integer, dst::Integer)
    """ This is the sum of densities between the anihilation and creation operators
        (including both endpoints) which gives the entire fermionic sign factor.
        The contributions from up and down states simply add.
    """
    if  src > dst
        return (-1)^count_occupancies(s, dst, src)
    end
    return (-1)^count_occupancies(s, src, dst)
end
