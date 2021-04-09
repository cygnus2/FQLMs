#===============================================================================

    hamiltonian_construction.jl - LR, January 2021

    Everything that is necessary for the construction of the Hamiltonian.

===============================================================================#
include("typedefs.jl")

function set_bits(ind::Array{LinkIndex,1})::LinkType
    """ Returns a state where the indicies are set on the specified links.
    """
    state = LinkType(0)
    for i in ind
        state = create(state, i)
    end
    return state
end

function apply_u(state::LinkType, p::Plaquette; sign::Bool=true)::Tuple{Union{Nothing,LinkType},Integer}
    """ Wrapper to apply U
    """
    return _apply_plaquette_operator(state, p, [false,false,true,true])
end

function apply_u_dagger(state::LinkType, p::Plaquette; sign::Bool=true)::Tuple{Union{Nothing,LinkType},Integer}
    """ Wrapper to apply U+
    """
    return _apply_plaquette_operator(state, p, [true,true,false,false])
end

function _apply_plaquette_operator(state::LinkType, plaquette::Plaquette, mask::Vector{Bool}; sign::Bool=true)::Tuple{Union{Nothing,LinkType},Integer}
    """ Applies the U operator

            c1+ c2+ c3 c4

        or the U+ term

            c1 c2 c3+ c4+

        to a given plaquette in a given state (link configuration starting
        with x-mu link and going counter-clockwise).

        Notes:
         -  The orientation (clockwise or anticlockwise) does not matter, the
            sign is in fact the same.
         -  For the summing of the occupancies, which gives the sign: in the
            Python version we first counted, then we changed the state. These
            operations commute though, because we always add the created/annihilated
            particles to the sum. Those are an even factor so the sign is preserved.
    """
    n = 0
    max_link = maximum(plaquette)
    new_state = state
    for (k,site) in enumerate(plaquette)
        new_state = mask[k] ? annihilate(new_state, site) : create(new_state, site)
        if !isnothing(new_state) && sign
            n += count_occupancies(new_state, plaquette[k], max_link)
        else
            return nothing, 0
        end
    end
    return new_state, (-1)^n
end


function do_single_state(state::LinkType, plaquettes::Array{Plaquette,1})::Array{Tuple{LinkType,LinkType,Integer},1}
    """ Constructs matrix-elements for a single Fock state, suitable for parallel
        computing.
    """
    states = Array{Tuple{LinkType,LinkType,Integer},1}()
    for p in plaquettes
        (new_state, s) = apply_u(state, p)

        # If U term was not successful, try the U^dagger term.
        # (the order could have been switched - there's always only one
        # possibility for overlap to be generated)
        if isnothing(new_state)
            (new_state, s) = apply_u_dagger(state, p)
        end

        # If a state was constructed, we'll add it to the list.
        if !isnothing(new_state)
            push!(states, (state, new_state, s))
        end
    end
    return states
end


function construct_hamiltonian(
    lookup_table::LookupDict,
    ilookup_table::InvLookupDict,
    latt::Lattice;
    silent::Bool=true
)
    """ Actually builds the Hamiltonian and returns a Hamiltonian in sparse
        representation.
    """
    # Loop through all Fock states and create the overlap matrix. First step:
    # do it naively (with some doubled work). Then try to improve on that (by
    # using, e.g., Hermiticity).
    row = Vector{HilbertIndex}()
    col = Vector{HilbertIndex}()
    data = Vector{DType}()

    # Get plaquettes we need to cycle.
    plaquettes = get_plaquettes(latt)

    # This counts the flippable plaquettes, which is essentially the construction
    # of the diagonal part of the Hamiltonian. This may be useful for later purposes
    # when opting for the computation of the fidelity susceptibility and related
    # quantities.
    n_flip_vector::Array{IType,1} = zeros(IType, length(lookup_table))

    # Loop through the list.
    for k=1:length(lookup_table)
        result = do_single_state(LinkType(lookup_table[k]), plaquettes)
        for (_, new_state, s) in result
            i = get(ilookup_table, new_state, nothing)
            if !isnothing(i)
                push!(row, k)
                push!(col, i)
                push!(data, s)
                n_flip_vector[k] += 1
            end
        end
    end
    return GaussLatticeHamiltonian(row, col, data, length(lookup_table)), n_flip_vector
end
