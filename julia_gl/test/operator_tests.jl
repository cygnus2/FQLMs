#===============================================================================

    operator_tests.jl - LR, January 2021

    Unittests for application of operators to the lattice. All operators are
    checked versus the Python implementation which is should be correct.

===============================================================================#
include("../src/typedefs.jl")
include("../src/operators/gl_operators.jl")
using Test

# Tests only for small systems.
const LinkType = SmallLinkState

@testset "julia operator tests" begin

    @testset "link state indexing" begin
        """ Ensures correct indexing of the states.
        """
        state = SmallLinkState(5)
        @test state[UInt8.([1,3,2])] == 3

        state = LargeLinkState(5)
        @test state[UInt8.([1,3,2])] == 3
    end

    @testset "parity checks 2x2" begin
        """ Make sure the parity operator is doing the correct thing for a
            simple 2x2 lattice.
        """
        # Make parity operator by Hand.
        op = BitmapOperator("parity", [3, 6, 1, 8, 7, 2, 5, 4])

        # Particles on
        state = SmallLinkState(5)
        new_state, sign = apply_operator(op, state)
        @test new_state == state && sign==-1
    end


    @testset "parity checks 2x2x4" begin
        """ Make sure the parity operator is doing the correct thing for a
            2x2x4 lattice.
        """
        # Make parity operator by Hand.
        op = gl_operators["2x2x4"]["parity"]

        # Exchange of two particles, parity eigenstate.
        state = SmallLinkState((1<<2)+(1<<38))
        new_state, sign = apply_operator(op, state)
        @test new_state == state && sign == -1

        # Two random particles, not an eigenstate.
        # [3, 1] --> [39, 4] (positive sign)
        initial = SmallLinkState((1<<2)+(1<<0))
        cstate, csign = apply_operator(op, initial)

        expected = SmallLinkState((1<<38)+(1<<3))
        @test cstate == expected && csign == 1
    end
end
