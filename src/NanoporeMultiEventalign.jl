__precompile__(true)

module NanoporeMultiEventalign
using BioSequences
using FASTX


include("utils.jl")
include("persistence.jl")
include("dtw.jl")



export dtw
export loadfasta

end # module
