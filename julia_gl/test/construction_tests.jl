#===============================================================================

    tests_2d.jl - LR, November 2020

    Some unittests for the HO stuff.

===============================================================================#
include("../src/typedefs.jl")
include("../src/hamiltonian_construction.jl")
using Test

@testset "hamiltonian construction tests" begin
    
    @testset "plaquette U operator" begin
        """ Check if the plaquette operator U does what it should do.
        """
        # Testcase 1: should give an overlap.
        # |0000 0001 0010 0000> ---> |0000000010010000>
        latt = LinkLattice([2,4])
        plaquettes = get_plaquettes(latt)

        state = parse(LinkState, "0000000100100000"; base=2)
        expected_state = parse(LinkState, "0000000010010000"; base=2)
        println(Int.(plaquettes[3].links))
        constructed_state, _ = apply_u(state, plaquettes[3])
        println(bitstring(state))
        println(bitstring(constructed_state))
        println(bitstring(expected_state))
        @test expected_state == constructed_state

        # # Testcase 2: should anihilate
        # # |0010 0000 0010 0000> ---> 0
        # state = int("0010 0000 0000 0000".replace(" ", ""), 2)
        # constructed_state, _ = apply_u(state, builder.plaquettes[2])
        # assert constructed_state == 0
        #
        # # ---
        #
        # # Testcase 3: 3D builder.
        # param["L"] = [2,2,2]
        # builder = HamiltonianBuilder(param, states=[])
        # state = set_bits([20,7])
        # expected = set_bits([19,14])
        #
        # p_ind = [19, 14, 7, 20]
        # constructed, _ = apply_u(state, p_ind + [set_bits(p_ind)])
        # assert expected == constructed
        #
        # constructed, _ = apply_u(state, builder.plaquettes[19])
        # assert expected == constructed
    end


end
