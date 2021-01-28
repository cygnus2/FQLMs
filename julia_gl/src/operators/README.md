# Operators
This code provides some functionality for the efficient application of operators to arbitrary states. The implementation rests on the mapping of bits from the initial to the final state. These maps could in principle be found via expressing the operator in question in second quantized form, however, this is challenging (I think) for rotation operators around arbitrary axes. Therefore, these maps are found with a suitable Python implementation and stored in in the `operator_masks.yml` file, from where they can be easily read. They are specific to the lattice size / geometry. 

## Bit maps
TBW
