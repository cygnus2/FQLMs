# Operators
This code provides some functionality for the efficient application of operators to arbitrary states. The implementation rests on the mapping of bits from the initial to the final state. These maps could in principle be found via expressing the operator in question in second quantized form, however, this is challenging (I think) for rotation operators around arbitrary axes. Therefore, these maps are found with a suitable Python implementation and exported such that they prodive a valid Julia dictionary. This is stored, in appropriate form, in `gl_operators.jl` which is the only thing we need to include to make use of them.

## Bit maps
TBW
