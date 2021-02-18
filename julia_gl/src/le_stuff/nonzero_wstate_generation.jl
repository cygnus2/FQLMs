#===============================================================================

    nonzero_wstate_generation.jl - LR, February 2021

    This is a script to create non-zero winding states from maximally flippable
    states in the [0,0,0] winding sector. The idea is to only alter these states
    slightly to introduce as litle deviation from the maximally flippable status
    as possible. One way to do this (I'm sure there are others) is to only
    change flip the spins along a straight line between the two sides of the box
    when the line is fully occupied or empty (constant occupancy/spin). For
    periodic boundary conditions, this is essentially a line/loop defect and the
    Gauss law is still intact. It turns out that this produces maximally
    flippable states in the winding sector with one flux along the direction
    of the line.

===============================================================================#

function find_lines(latt::LinkLattice, dir::Direction)::Set{Array{LinkType,1},1}
    """ Finds the straight lines that connect two sides of the lattice in a given
        direction. These are represented by a list of link indicies.

        For instance on a 2x2x2 lattice we may follow the y-direction starting
        from the lower left corner. The corresponding links would be [2, 8]. Note
        that those are 1-based (Julia default).
    """
    verts = get_vertices(latt)

    # Lazy version: simply loop over all vertices, find the lines in the specified
    # direction and then remove the doubles. This way we don't have to find the
    # vertices in a given plane. Doubles will not be stored because we use a set.
    lines = Set{Array{Linktype,1},1}()
    for i=1:length(verts)
        v0 = verts[i] # We need the initial vertex.
        next = verts[i]

        line = Array{LinkType,1}()
        while next[]
            push(line, )

        end


        _shift_index()

        # Add the sorted line (to avoid doubles) to the set.
        push!(lines, sort(line))
    end
    return lines
end


function increase_winding(base_states::Array{LinkType,1}, latt::Lattice, dir::Direction; increment::IType=1)::Array{LinkType,1}
    """ Takes in a list of states and returns states with winding changed by one
        unit (positive or negative, specified by increment) along a given
        direction (x=1, y=2, z=3).
    """
    # The lines along which flipping spins is attempted.
    lines = find_lines(latt, dir)

    # Using a set immediately lets us prevent doubles (on the cheap!)
    winding_states = Set{LinkType}()
    for state in base_states
        for line in lines
            new_state = state

            # Try to flip all links in a given line. As soon one doesn't work,
            # the line is not flippable.
            for link in line
                new_state = increment > 0 ? create(new_state,link) : annihilate(new_state,link)
                if isnothing(new_state)
                    break
                end
            end
            if !isnothing(new_state)
                push!(winding_states, new_state)
            end
        end
    end

    # We want an array, not a set.
    return sort(collect(winding_states))
end
