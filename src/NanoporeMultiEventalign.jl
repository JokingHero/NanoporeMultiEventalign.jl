__precompile__(true)

module NanoporeMultiEventalign
using BioSequences
using FASTX
using Plots

include("utils.jl")
include("persistence.jl")
include("dtw.jl")


export dtw
export loadfasta, loadkmers, fasta_to_kmer_valus

end # module
