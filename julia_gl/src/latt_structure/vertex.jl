#===============================================================================

    vertex.jl - LR, Jan 2021

    Representation of a vertex is done via a plain array.

    Contents of a vertex:
            1  2  3  4  5  6 ...
            x -x  y -y  z -z ...

    The representation could be done in any other way (notably, with a struct),
    this file is there to accommodate potential future implementations.

===============================================================================#
# const Vertex = Vector{LinkIndex}

struct Vertex
    links::Vector{LinkIndex} # Lattice dimensions [Lx, Ly, (Lz)].
    i::Union{SiteIndex,Nothing} # Index in the lattice. Not strictly necessary, but convenient.

    Vertex() = new([],nothing)
    Vertex(links::Vector{Integer}) = new(LinkIndex.(links),nothing)
    Vertex(links::Vector{T}, i::Integer) where T<:Integer = new(LinkIndex.(links),SiteIndex(i))
end
Base.getindex(v::Vertex, i::Integer)::LinkIndex = v.links[i]
Base.push!(v::Vertex, l::Integer) = push!(v.links, LinkIndex(l))


Base.firstindex(v::Vertex) =  1
Base.lastindex(v::Vertex) = length(v.links)

# For iteration.
Base.length(v::Vertex) = length(v.links)
function Base.iterate(v::Vertex, state::Integer=1)
    state > length(v) ? nothing : (v.links[state], state+1)
end

# For pretty printing.
Base.show(io::IO, v::Vertex) = print(io, "Vertex [(i=",v.i,"), links=", Int.(v), "]")
