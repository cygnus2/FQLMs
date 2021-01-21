#===============================================================================

    lattice.jl - LR, January 2021

    Type that represents a lattice - essentially gives the code the underlying
    structure to operate on.

===============================================================================#
using IterTools


struct LinkLattice <: Lattice
    L::Vector{IType} # Lattice dimensions [Lx, Ly, (Lz)].
    S::Vector{IType} # Shift operators to get the link indices.
    d::Integer # Spatial dimension of the lattice

    LinkLattice(L::Vector{IType}) = begin
        S::Array{IType,1} = [1]
        for l in L
            push!(S, S[end]*l)
        end
        new(L,S,length(L))
    end
end


function _shift_index(i::SiteIndex, dir::Integer, l::LinkLattice)::SiteIndex
    """ Shifts an grid index into direction d under consideration of  PBC.
        Works only for positive shifts (in positive coordinate axes) but
        this is also the only use case.

        Notes:
            -   dir is 1-based (1=x, 2=y, 3=z, ... )
            -   Index is now 1-based, howver, logic still from zero-based
                indicies (Python) but adapted
    """
    n = i-1 + l.S[dir]
    if mod(n,l.S[dir+1]) < l.S[dir]
        n -= l.S[dir+1] #- S[d-1]
    end
    return n+1
end


function _get_single_vertex(latt::LinkLattice, i::SiteIndex)::Vertex
    """ Returns the indices [+x, -x, +y, -y, ...] of the links for the i-th
        vertex *on the full lattice* (not the sublattice) under consideration
        of periodic boundary conditions.

        Note: works only for L^d lattices for now.
    """
    # Shorthand to avoid self all the time.
    # Index in bit string.
    j = latt.d*(i-1) + 1

    # Loop through the dimensions and add the indicies of + and - directions.
    vert = Vertex()
    for k = 1:latt.d
        # Step forward is always 'on site'.
        push!(vert, j+k-1)

        # Sep backward in k direction.
        l = j + k - 1 - latt.d*latt.S[k]
        if mod(i-1,latt.S[k+1]) < latt.S[k]
            l = l + latt.d*latt.S[k+1]
        end
        push!(vert, l)
    end
    return vert
end


function get_vertices(latt::LinkLattice)::Array{Vertex,1}
    """ Returns a list of vertices.
    """
    return [_get_single_vertex(latt,k) for k=SiteIndex.(1:latt.S[end])]
end


function get_plaquettes(latt::LinkLattice; separate_lists::Bool=false)
    """ Gets the entire list of plaquettes for a given lattice.
    """
    # Get a list of vertices to work with.
    v = get_vertices(latt)

    # Find plaquettes by looping over all grid points.
    (p_xy, p_yz, p_xz) = Array{Plaquette,1}(), Array{Plaquette,1}(), Array{Plaquette,1}()
    for n = SiteIndex.(1:latt.S[end])
        # xy plane.
        j = _shift_index(_shift_index(n, 1, latt), 2, latt) # shifted by Sx and Sy
        push!(p_xy, Plaquette([v[n][1], v[j][4], v[j][2], v[n][3]]))

        # In 3D, we have two additional plaquettes.
        if latt.d == 3
            # yz plane
            j = _shift_index(_shift_index(n, 3, latt), 2, latt) # shifted by Sy and Sz
            push!(p_yz, Plaquette([v[n][3], v[j][6], v[j][4], v[n][5]]))

            # xz plane.
            j = _shift_index(_shift_index(n, 1, latt), 3, latt) # shifted by Sx and Sz
            push!(p_xz, Plaquette([v[n][1], v[j][6], v[j][2], v[n][5]]))
        end
    end

    # Check if the right amount of plaquettes was found and if so, return
    # the list.
    plaquettes = vcat(p_xy, p_yz, p_xz)
    @assert length(plaquettes) == (2^(latt.d-1) -1) * latt.S[end]

    if separate_lists
        return p_xy, p_yz, p_xz
    end
    return plaquettes
end


function winding_sectors(latt::LinkLattice)
    """ A generator for all winding number sectors.
    """
    winding_numbers = nothing
    if latt.d == 2
        winding_numbers = reverse(latt.L).+1
    elseif latt.d == 3
        winding_numbers = (
            latt.L[2]*latt.L[3] + 1,
            latt.L[1]*latt.L[3] + 1,
            latt.L[1]*latt.L[2] + 1
        )
    else
        error("Only 2D and 3D lattices are allowed.")
    end

    tups = product(map(n->(0:(n-1)), winding_numbers)...)
    return collect.(reshape(collect(tups), (length(tups),1)))
end
