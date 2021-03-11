#===============================================================================

    tests_2d.jl - LR, November 2020

    Some unittests for the HO stuff.

===============================================================================#
include("../src/typedefs.jl")
const LinkType = SmallLinkState
include("../src/hamiltonian_construction.jl")
using Test

# All tests here are done with the smaller representation.


@testset "hamiltonian construction tests" begin

    @testset "plaquette U operator" begin
        """ Check if the plaquette operator U does what it should do.
        """
        # Testcase 1: should give an overlap.
        # |0000 0001 0010 0000> ---> |0000000010010000>
        latt = LinkLattice([2,4])
        plaquettes = get_plaquettes(latt)

        state = parse(LinkType, "0000000100100000"; base=2)
        expected_state = parse(LinkType, "0000000010010000"; base=2)
        constructed_state, _ = apply_u(state, plaquettes[3])
        # println(Int.(plaquettes[3].links))
        # println(bitstring(state))
        # println(bitstring(constructed_state))
        # println(bitstring(expected_state))
        @test expected_state == constructed_state

        #---
        # Testcase 2: should anihilate
        # |0010 0000 0010 0000> ---> 0
        state = parse(LinkType, replace("0010 0000 0000 0000", " "=>""); base=2)
        constructed_state, _ = apply_u(state, plaquettes[2])
        @test isnothing(constructed_state)

        constructed_state, _ = apply_u_dagger(state, plaquettes[4])
        @test isnothing(constructed_state)

        # ---
        # Testcase 3: 2x2x2
        latt = LinkLattice([2,2,2])

        state = set_bits(LinkIndex.([21, 8]))
        expected = set_bits(LinkIndex.([20, 15]))
        p = Plaquette(LinkIndex.([20, 15, 8, 21]))
        #
        constructed, _ = apply_u(state, p)
        @test expected == constructed
        #
        constructed, _ = apply_u(state, get_plaquettes(latt)[15])
        @test expected == constructed

        # ---
        # Testcase 4: 2x2x4
        latt = LinkLattice([2,2,4])
        plaquettes = get_plaquettes(latt)

        # This is an unflippable state, so it should annihilate for all operators.
        state = LinkType(11420033508079)
        for p in plaquettes
            constructed, _ = apply_u(state, p)
            @assert isnothing(constructed)
            constructed, _ = apply_u_dagger(state, p)
            @assert isnothing(constructed)
        end
        @test true

        # Flippable state, fifth plaquette should be flippable by U.
        state = parse(LinkType, "100101100101100001001000101100110100010110110001"; base=2)
        expected = parse(LinkType, "100101100101100001001000100100101110010110110001"; base=2)
        constructed, _ = apply_u(state, plaquettes[6])
        @test expected == constructed
    end


    @testset "plaquette signs" begin
        """ Check signs of plaquettes. This is tested on a maximally flippable
            state on a 2x2x2 lattice and a few random plaquettes.

            As opposed to above, this includes the signs.
        """
        latt = LinkLattice([2,2,2])
        plaquettes = get_plaquettes(latt)
        state = LinkState([22, 21, 20, 18, 14, 13, 12, 11, 7, 5, 4, 3])


        # First plaquette.
        # p1 = Plaquette(LinkIndex.([2, 9, 14, 3]))
        p1 = plaquettes[9]
        e_state = LinkState([22, 21, 20, 18, 13, 12, 11, 9, 7, 5, 4, 2])
        e_sign = -1
        c_state, c_sign = apply_u(state, p1)
        @test c_state == e_state && c_sign == e_sign

        # Second plaquette.
        # p2 = Plaquette(LinkIndex.([19, 23, 13, 20]))
        p2 = plaquettes[7]
        e_state = LinkState([23, 22, 21, 19, 18, 14, 12, 11, 7, 5, 4, 3])
        e_sign = -1
        c_state, c_sign = apply_u(state, p2)
        @test c_state == e_state && c_sign == e_sign

        # Third plaquette.
        # p3 = Plaquette(LinkIndex.([13, 18, 1, 15]))
        p3 = plaquettes[21]
        e_state = LinkState([22, 21, 20, 15, 14, 12, 11, 7, 5, 4, 3, 1])
        e_sign = -1
        c_state, c_sign = apply_u_dagger(state, p3)
        @test c_state == e_state && c_sign == e_sign
    end
end
