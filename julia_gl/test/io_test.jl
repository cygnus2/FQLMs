#===============================================================================

    io_test.jl - LR, November 2020

===============================================================================#
include("../src/typedefs.jl")
include("../src/io/data_import.jl")
using Test

@testset "julia import tests" begin

    @testset "long type conversion" begin
        """ HDF5 does not store larger objets than Int64, therefore we have to
            store larger numbers as two integers. This routines tests the conversion.
        """
        original = LargeLinkState.([1074260571646206819420, 1089901011134353970538,  1526095490192680073940, 1598215294499136381873, 2123163061129091160270, 2179642425947400317085, 2542724056922244896610, 2599203421740554053425, 3124151188370508831822, 3196270992676965139755, 3632465471735291243157, 3648105911223438394275])

        # Make back-and-forth conversions, those should be the same.
        converted = _convert_links_to_HDF5(original)
        recreated = _convert_links_from_HDF5(converted)
        @test all(recreated.==original)
    end
end
