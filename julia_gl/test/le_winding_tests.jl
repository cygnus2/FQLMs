#===============================================================================

    le_winding_tests.jl - LR, February 2021

    Some unittests for the HO stuff.

===============================================================================#
include("../src/typedefs.jl")
const LinkType = SmallLinkState

include("../src/le_stuff/nonzero_wstate_generation.jl")
using Test

@testset "julia le tests" begin

    @testset "find lines" begin
        """ Make sure the correct lines along the directions are found.
        """
        latt = LinkLattice([2,2,2])

        exp_xlines = Set([[1,4], [7,10], [13, 16], [19,22]])
        cons_xlines = find_lines(latt, xpos)
        @test length(setdiff(exp_xlines, cons_xlines)) == 0

        exp_ylines = Set([[2,8], [5,11], [14,20], [17,23]])
        cons_ylines = find_lines(latt, ypos)
        @test length(setdiff(exp_ylines, cons_ylines)) == 0

        exp_zlines = Set([[3,15], [6,18], [9,21], [12,24]])
        cons_zlines = find_lines(latt, zpos)
        @test length(setdiff(exp_zlines, cons_zlines)) == 0
    end

    @testset "find maximally flippable winding states" begin
        """ Make sure we find the correct windin states.
        """
        latt = LinkLattice([2,2,2])

        # These are the maximally flippable states for 2x2x2.
        base_states_222 = [3816540, 3872106, 5421780, 5678001, 7542990, 7743645,
                            9033570, 9234225, 11099214, 11355435, 12905109, 12960675]


        # Expected [1,0,0] states for a 2x2x2 lattice.
        expected_1_0_0 = [5421789, 7744221, 7780509, 7781076, 9033579, 11356011, 11392299, 11392866]
        constructed_1_0_0 = increase_winding(LinkType.(base_states_222), latt, xpos; increment=1)
        for k=1:length(expected_1_0_0)
            @assert expected_1_0_0[k] == constructed_1_0_0[k]
        end
        @test true


        # Expected [-1,0,0] states for a 2x2x2 lattice.
        expected_m1_0_0 = [5384349, 5384916, 5421204, 7743636, 8996139, 8996706, 9032994, 11355426]
        constructed_m1_0_0 = increase_winding(LinkType.(base_states_222), latt, xpos; increment=-1)
        for k=1:length(expected_m1_0_0)
            @assert expected_m1_0_0[k] == constructed_m1_0_0[k]
        end
        @test true
    end
end
