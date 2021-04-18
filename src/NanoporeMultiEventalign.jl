__precompile__(true)

module NanoporeMultiEventalign
using BioSequences
using FASTX
using HDF5
using Changepoints
using DynamicAxisWarping
using ChangePointDetection

include("utils.jl")
include("persistence.jl")
include("dtw.jl")

X = [1.3,3.3,4.3,5.3]
C = [3.3,04.3,5.3,6.2]

println(lsdd(X,C))

export nanopore_dtw
export loadnanoporefast5, loadfasta, loadkmers, fasta_to_kmer_values
export kmerDist, Bhattacharyya


end # module
