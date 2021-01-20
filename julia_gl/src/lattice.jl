#===============================================================================

    lattice.jl - LR, January 2021

    Type that represents a lattice - essentially gives the code the underlying
    structure to operate on.

===============================================================================#
abstract type Lattice end


struct LinkLattice <: Lattice
    L::Vector{Int} # Lattice dimensions [Lx, Ly, (Lz)].
    S::Vector{Int} # Shift operators to get the link indices.
    d::Integer # Spatial dimension of the lattice

    LinkLattice(L::Vector{Int}) = begin
        S::Array{Int,1} = [1]
        for l in L
            push!(S, S[end]*l)
        end
        new(L,S,length(d))
    end
end


function get_vertex_links(latt::LinkLattice, i::Integer)
    """ Returns the indices [+x, -x, +y, -y, ...] of the links for the i-th
        vertex *on the full lattice* (not the sublattice) under consideration
        of periodic boundary conditions.

        Note: works only for L^d lattices for now.
    """
    # Shorthand to avoid self all the time.
    # Index in bit string.
    j = latt.d*i

    # Loop through the dimensions and add the indicies of + and - directions.
    ind = []
    for k = 1:latt.d
        # Step forward is always 'on site'.
        push!(vert, j+k)

        # Sep backward in k direction.
        l = j + k - latt.d*lat..S[k]
        if mod(i,latt.S[k]) < latt.S[k]
            l = l + latt.d*latt.S[k+1]
        end
        push!(vert, l)
    end
    return vert
end

struct get_plaquettes(latt::LinkLattice)
    """ Produces a list of plaquettes that represent the lattice. This is a
        list of plaquettes represented as

            plaquette = [link_1, link_2, link_3, link_4]

        where link link_1 and link_2 as well as link_3 and link_4 are
        connected through the hopping like

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
    """
    # Find plaquettes by looping over all grid points.
    p_xy, p_yz, p_xz = Array{Plaquette,1}(), Array{Plaquette,1}(), Array{Plaquette,1}()

    for n in range(latt.S[end])
        # Contents of a vertex:
        # 0  1  2  3  4  5 ...
        # x -x  y -y  z -z ...
        vn = get_vertex_links(latt, n)

        # xy plane.
        j = self.shift_index(self.shift_index(n, 0), 1) # shifted by Sx and Sy
        vn_xy = self.get_vertex_links(j)
        ind = [vn[0], vn_xy[3], vn_xy[1], vn[2]]
        push!(p_xy, Plaquette(ind))

        # In 3D, we have two additional plaquettes.
        if self.d == 3
            # yz plane
            j = self.shift_index(self.shift_index(n, 2), 1) # shifted by Sy and Sz
            vn_yz = self.get_vertex_links(j)
            ind=[vn[2], vn_yz[5], vn_yz[3], vn[4]]
            p_yz.append(ind + [set_bits(ind)])

            # xz plane.
            j = self.shift_index(self.shift_index(n, 0), 2) # shifted by Sx and Sz
            vn_xz = self.get_vertex_links(j)
            ind = [vn[0], vn_xz[5], vn_xz[1], vn[4]]
            p_xz.append(ind + [set_bits(ind)])
        end
    end

    # Check if the right amount of plaquettes was found and if so, return
    # the list.
    plaquettes = p_xy + p_yz + p_xz
    assert len(plaquettes) == (2**(self.d-1) -1) * S[-1]

    if separate_lists
        return p_xy, p_yz, p_xz
    end
    return plaquettes
end
