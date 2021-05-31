__precompile__(true)

module NanoporeMultiEventalign
using BioSequences
using FASTX
using HDF5
using Plots

include("utils.jl")
include("persistence.jl")
include("dtw.jl")


export nanopore_dtw, multi_nanopore_dtw, nanopore_dtw_plot
export loadnanoporefast5, loadfasta, loadkmers, fasta_to_kmer_values
export kmerdist, bhattacharyya, multi_bhattacharyya

end # module
