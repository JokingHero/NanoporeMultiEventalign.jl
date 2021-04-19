__precompile__(true)

module NanoporeMultiEventalign
using BioSequences
using FASTX
using HDF5

include("utils.jl")
include("persistence.jl")
include("dtw.jl")



export nanopore_dtw
export loadnanoporefast5, loadfasta, loadkmers, fasta_to_kmer_values
export kmerdist, bhattacharyya


end # module
