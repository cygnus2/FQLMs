#===============================================================================

    state.jl - LR, Jan 2021

    Representation of states is done via an Integer.

===============================================================================#
# That's the occupation state of a lattice of links.
const LinkIndex = UInt8
const SmallLinkState = UInt64
const LargeLinkState = UInt128
const LinkState = Union{SmallLinkState,LargeLinkState}


function create(s::LinkState, i::Integer)::Union{Nothing,LinkState}
    """ Destroys a particle at position i, false if there's a particle already.
    """
    # Typeof here is *very* important, otherwise the 1 defaults to Int64 and it
    # won't be possible to use larger representations than 63 bits.
    shift = typeof(s)(1) << (i-1)
    if s & shift == 0
        return s ⊻ shift
    else
        return nothing
    end
end

function annihilate(s::LinkState, i::Integer)::Union{Nothing,LinkState}
    """ Destroys a particle at position i, false if there's no particle to destroy.
    """
    shift = typeof(s)(1) << (i-1)
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
    for i=src-1:dst-1
        c += ((s >>> i) & 1)
    end
    return c
end


function sign(s::LinkState, src::Integer, dst::Integer)
    """ This is the sum of densities between the anihilation and creation operators
        (including both endpoints) which gives the entire fermionic sign factor.
        The contributions from up and down states simply add.
    """
    println("WHAT? This shouldn't be here..")
    if  src > dst
        return (-1)^count_occupancies(s, dst, src)
    end
    return (-1)^count_occupancies(s, src, dst)
end


# Indexing.
function Base.getindex(s::LinkState, i::Unsigned)::LinkState
    """ Note: Indexing works with the lowest bits first (i.e., the ones from the
        right). I guess it doesn't matter how it is used, only the implementation
        here needs to be consistent.

        Only positive integers are allowed. Negatives would only give zeros.
    """
    return s >>> (i-1) & 1
end
Base.firstindex(p::LinkState) =  1


function Base.getindex(s::T, iarr::Array{S,1})::T where T<:LinkState where S<:Unsigned
    """ Array indexing shuffles the bits around according to the indicies.
        Example: Here we want to flip the first two bits:

            (001101)[2,1,3,4,5,6] = (001110)

        Can this be improved somehow?
    """
    perm::T = 0
    for (k, i) in enumerate(iarr)
        perm += s[i] << (k-1)
    end
    return perm
end
