#===============================================================================

    le_winding_tests.jl - LR, February 2021

    Some unittests for the HO stuff.

===============================================================================#
include("../src/typedefs.jl")
include("../src/le_stuff/nonzero_wstate_generation.jl")
using Test

@testset "julia le tests" begin

    @testset "find lines" begin
        """ Make sure the correct lines along the directions are found.
        """
        latt = LinkLattice([2,2,2])

        exp_xlines = Set([[1,2], [7,10], [13, 16], [19,22]])
        cons_xlines = find_lines(latt, xpos)
        @test length(setdiff(exp_xlines, cons_xlines)) == 0

        exp_ylines = Set([[2,8], [5,11], [14,20], [17,23]])
        cons_ylines = find_lines(latt, ypos)
        @test length(setdiff(exp_ylines, cons_ylines)) == 0

        exp_zlines = Set([[3,15], [6,18], [9,21], [12,24]])
        cons_zlines = find_lines(latt, zpos)
        @test length(setdiff(exp_zlines, cons_zlines)) == 0
    end
end
