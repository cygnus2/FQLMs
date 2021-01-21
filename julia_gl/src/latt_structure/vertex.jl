#===============================================================================

    vertex.jl - LR, Jan 2021

    Representation of a vertex is done via a plain array.

    Contents of a vertex:
            1  2  3  4  5  6 ...
            x -x  y -y  z -z ...

    The representation could be done in any other way (notably, with a struct),
    this file is there to accommodate potential future implementations.

===============================================================================#
const Vertex = Vector{LinkIndex}
