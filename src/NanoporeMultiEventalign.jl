__precompile__(true)

module NanoporeMultiEventalign
using BioSequences
using FASTX
using Plots

include("persistence.jl")
include("dtw.jl")
include("utils.jl")


export dtw
export loadfasta, loadkmers, plotdata

end # module
