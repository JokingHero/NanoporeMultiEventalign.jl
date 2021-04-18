__precompile__(true)

module NanoporeMultiEventalign
using BioSequences
using FASTX
using HDF5
using DynamicAxisWarping
using ChangePointDetection

include("utils.jl")
include("persistence.jl")
include("dtw.jl")
include("ChangePointDetection.jl")

ts1 = rand(50)
ts2 = 1.5 * rand(50) .+ 1
ts = vcat(ts1, ts2)
profile = lsdd_profile(ts; window = 1)
points = getpoints(profile, threshold = 1.3)
print(points)

export nanopore_dtw
export loadnanoporefast5, loadfasta, loadkmers, fasta_to_kmer_values
export kmerDist, Bhattacharyya


end # module
