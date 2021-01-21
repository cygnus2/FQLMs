#===============================================================================

    tests_2d.jl - LR, November 2020

    Some unittests for the HO stuff.

===============================================================================#
include("../src/typedefs.jl")
using Test

@testset "julia gauss lattice tests" begin

    @testset "vertex creation" begin
        # Make data.
        latt = LinkLattice([2,2,4])
        vertices = get_vertices(latt)

        # All vertices to check.
        sites = [1 5 10 16]
        expected = [
            [1, 4, 2, 8, 3, 39],
            [13, 16, 14, 20, 15, 3],
            [28, 25, 29, 35, 30, 18],
            [46, 43, 47, 41, 48, 36]
        ]
        for (n, site) in enumerate(SiteIndex.(sites))
            # println("--- $site ---")
            for dir = 1:2*latt.d
                # println(expected[n][dir] , " // ", vertices[site][dir])
                @assert expected[n][dir] == vertices[site][dir]
            end
        end
        @test true
    end

    @testset "PBC index shift 3D" begin
        """ Check if PBC are enforced correctly.
        """
        expected = [
            [2, 3, 5],
            [1, 4, 6],
            [4, 1, 7],
            [3, 2, 8],

            [6, 7, 9],
            [5, 8, 10],
            [8, 5, 11],
            [7, 6, 12],

            [10, 11, 13],
            [9, 12, 14],
            [12, 9, 15],
            [11, 10, 16],

            [14, 15, 1],
            [13, 16, 2],
            [16, 13, 3],
            [15, 14, 4]
        ]

        latt = LinkLattice([2,2,4])
        for n = SiteIndex.(1:latt.S[end])
            for dir = 1:latt.d
                @assert expected[n][dir] == _shift_index(n,dir,latt)
            end
        end
        @test true
    end


    @testset "plaquettes 3D" begin
        """ Check if the correct list of plaquettes is found in 3D.
        """
        latt = LinkLattice([2,2,2])
        plaquettes = get_plaquettes(latt)

        # Set up test case.
        # (these are zero-based indicies copied from Python, hence the +1 below)
        expected = [
            [0, 4, 6, 1],
            [3, 1, 9, 4],
            [6, 10, 0, 7],
            [9, 7, 3, 10],
            [12, 16, 18, 13],
            [15, 13, 21, 16],
            [18, 22, 12, 19],
            [21, 19, 15, 22],

            [1, 8, 13, 2],
            [4, 11, 16, 5],
            [7, 2, 19, 8],
            [10, 5, 22, 11],
            [13, 20, 1, 14],
            [16, 23, 4, 17],
            [19, 14, 7, 20],
            [22, 17, 10, 23],

            [0, 5, 12, 2],
            [3, 2, 15, 5],
            [6, 11, 18, 8],
            [9, 8, 21, 11],
            [12, 17, 0, 14],
            [15, 14, 3, 17],
            [18, 23, 6, 20],
            [21, 20, 9, 23],
        ]

        # Chekc the plaquette index.
        for (i, plaquette) in enumerate(plaquettes)
            # println("($i) ", Int.(plaquette), " // ", expected[i] .+ 1)
            for (j, link) in enumerate(plaquette)
                @assert expected[i][j] == link-1
            end
        end

        # # Check the mask array.
        # for i, p in enumerate(builder.plaquettes):
        #     mask = 0
        #     for k in expected[i]:
        #         mask = mask + (1 << k)
        #     assert mask == p[-1]

        @test true
    end



    @testset "plaquette U operator" begin
        """ Check if the plaquette operator U does what it should do.
        """
        # Testcase 1: should give an overlap.
        # |0000 0001 0010 0000> ---> |0000000010010000>
        latt_24 = LinkLattice([2,4])
        state = int('0000000100100000', 2)
        expected_state = int('0000000010010000', 2)
        constructed_state, _ = apply_u(state, builder.plaquettes[2])
        assert expected_state == constructed_state

        # Testcase 2: should anihilate
        # |0010 0000 0010 0000> ---> 0
        state = int('0010 0000 0000 0000'.replace(' ', ''), 2)
        constructed_state, _ = apply_u(state, builder.plaquettes[2])
        assert constructed_state == 0

        # ---

        # Testcase 3: 3D builder.
        param['L'] = [2,2,2]
        builder = HamiltonianBuilder(param, states=[])
        state = set_bits([20,7])
        expected = set_bits([19,14])

        p_ind = [19, 14, 7, 20]
        constructed, _ = apply_u(state, p_ind + [set_bits(p_ind)])
        assert expected == constructed

        constructed, _ = apply_u(state, builder.plaquettes[19])
        assert expected == constructed
    end


end
