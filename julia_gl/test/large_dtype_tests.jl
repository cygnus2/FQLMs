#===============================================================================

    tests_2d.jl - LR, November 2020

    Some unittests for the HO stuff.

===============================================================================#
include("../src/typedefs.jl")
include("../src/hamiltonian_construction.jl")
using Test

@testset "julia datatype tests" begin

    @testset "plaquette operator" begin
        """ Make sure, the bit representation can deal with large states.
        """
        latt = LinkLattice([2,2,6])

        state = LargeLinkState(52378704074343307122)
        p = Plaquette(LinkIndex.([56, 51, 68, 57]))

        new_state, sign = apply_u_dagger(state, p)
        @test new_state == LargeLinkState(199987559561131841394)

    end
end
