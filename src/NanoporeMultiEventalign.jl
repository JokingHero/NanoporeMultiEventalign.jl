__precompile__(true)

module NanoporeMultiEventalign
using HDF5

include("utils.jl")
include("persistence.jl")
include("dtw.jl")

export dtw

end # module
