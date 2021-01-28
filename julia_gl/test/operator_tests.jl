#===============================================================================

    operator_tests.jl - LR, January 2021

    Unittests for application of operators to the lattice. All operators are
    checked versus the Python implementation which is should be correct.

===============================================================================#
include("../src/typedefs.jl")
using Test

@testset "julia operator tests" begin

    @testset "link state indexing" begin
        """ Ensures correct indexing of the states.
        """
        state = SmallLinkState(5)
        @test state[UInt8.([1,3,2])] == 3

        state = LargeLinkState(5)
        @test state[UInt8.([1,3,2])] == 3
    end

    @testset "parity checks" begin
        """ Make sure the parity operator is doing the correct thing.
        """
        state = SmallLinkState(1232342)

        
    end
end
