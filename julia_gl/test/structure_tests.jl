#===============================================================================

    tests_2d.jl - LR, November 2020

    Some unittests for the HO stuff.

===============================================================================#
include("../src/typedefs.jl")
using Test

# All tests here are done with the smaller representation.
const LinkType = SmallLinkState

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


    @testset "plaquettes 2x2x2" begin
        """ Check if the correct list of plaquettes is found in a 2x2x2 lattice.
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


    @testset "plaquettes 2x2x4" begin
        """ Check if the correct list of plaquettes is found in a 2x2x4 lattice.
        """
        latt = LinkLattice([2,2,4])
        plaquettes = get_plaquettes(latt)

        # Set up test case.
        # (these are zero-based indicies copied from Python, hence the +1 below)
        expected = [
            [1, 5, 7, 2],
            [4, 2, 10, 5],
            [7, 11, 1, 8],
            [10, 8, 4, 11],
            [13, 17, 19, 14],
            [16, 14, 22, 17],
            [19, 23, 13, 20],
            [22, 20, 16, 23],
            [25, 29, 31, 26],
            [28, 26, 34, 29],
            [31, 35, 25, 32],
            [34, 32, 28, 35],
            [37, 41, 43, 38],
            [40, 38, 46, 41],
            [43, 47, 37, 44],
            [46, 44, 40, 47],

            [2, 9, 14, 3],
            [5, 12, 17, 6],
            [8, 3, 20, 9],
            [11, 6, 23, 12],
            [14, 21, 26, 15],
            [17, 24, 29, 18],
            [20, 15, 32, 21],
            [23, 18, 35, 24],
            [26, 33, 38, 27],
            [29, 36, 41, 30],
            [32, 27, 44, 33],
            [35, 30, 47, 36],
            [38, 45, 2, 39],
            [41, 48, 5, 42],
            [44, 39, 8, 45],
            [47, 42, 11, 48],

            [1, 6, 13, 3],
            [4, 3, 16, 6],
            [7, 12, 19, 9],
            [10, 9, 22, 12],
            [13, 18, 25, 15],
            [16, 15, 28, 18],
            [19, 24, 31, 21],
            [22, 21, 34, 24],
            [25, 30, 37, 27],
            [28, 27, 40, 30],
            [31, 36, 43, 33],
            [34, 33, 46, 36],
            [37, 42, 1, 39],
            [40, 39, 4, 42],
            [43, 48, 7, 45],
            [46, 45, 10, 48],
        ]

        # Chekc the plaquette index.
        for (i, exp) in enumerate(expected)
            # println("($i) ", Int.(plaquettes[i]), " // ", exp)
            for (j, link) in enumerate(exp)
                @assert plaquettes[i][j] == link
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
end
