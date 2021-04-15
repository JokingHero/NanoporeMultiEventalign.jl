__precompile__(true)

module NanoporeMultiEventalign
using BioSequences
using FASTX
using HDF5
using Changepoints
using DynamicAxisWarping

include("utils.jl")
include("persistence.jl")
include("dtw.jl")

a = [kmerDist(1,2),kmerDist(2,2),kmerDist(3,2),kmerDist(4,2),kmerDist(5,2)]
b = [kmerDist(3,2),kmerDist(4,2),kmerDist(5,2),kmerDist(6,2),kmerDist(7,2)]
println(dtw2(a,b))



export nanopore_dtw
export loadnanoporefast5, loadfasta, loadkmers, fasta_to_kmer_values
export kmerDist, Bhattacharyya


end # module
