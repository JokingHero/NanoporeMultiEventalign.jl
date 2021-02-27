__precompile__(true)

module NanoporeMultiEventalign
using BioSequences


include("utils.jl")
include("persistence.jl")
include("dtw.jl")
include("loadReferanceFa.jl")

export dtw

end # module
