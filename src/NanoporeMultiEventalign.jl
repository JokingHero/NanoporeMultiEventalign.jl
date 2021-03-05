__precompile__(true)

module NanoporeMultiEventalign

using HDF5
using BioSequences
using FASTX

include("utils.jl")
include("persistence.jl")
include("dtw.jl")

export dtw
export loadnanoporefast5, loadfasta

end # module
