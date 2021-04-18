
"""
Gets the corosponding mean and stdv from a fasta sequense using kmer data.
Takes in fasta data loaded from loadfasta and kmers loaded from loadkmers.
Returns a list of tuples where the values are {mean, stdv}.
"""
function fasta_to_kmer_values(fastadata, kmers)
    datastring = convert(String, fastadata[1])
    fastameans = []
    # gets the length of a single kmer key
    keylength = length(collect(keys(kmers))[1])
    for i in 1:(length(datastring)-keylength+1)
        kmerkey = datastring[i:(i+keylength-1)]
        #collects all the kmer means and stdv into a list of tuples
        push!(fastameans, tuple(parse(Float64, kmers[kmerkey][1]), parse(Float64, kmers[kmerkey][2])))
    end
    return fastameans
end

struct kmerDist
    mean::Float32
    sd::Float32
end


"""
computes the Bhattacharyya distance of two kmerDist x and y
"""
function Bhattacharyya(x::kmerDist,y::kmerDist)
    arg1 = 1/4*((x.sd*x.sd)/(y.sd*y.sd) + (y.sd*y.sd)/(x.sd*x.sd) + 2)
    arg2 = 1/4*((x.mean - y.mean)*(x.mean - y.mean)/(x.sd*x.sd + y.sd*y.sd))
    return 1/4*log(arg1 + arg2)
end


function dtw2(seq1::Array{kmerDist}, seq2::Array{kmerDist})
    D = dtw_cost_matrix2(seq1,seq2)
    return trackback2(D)
end



function pairwise3(s1::Array{kmerDist},s2::Array{kmerDist})
    [Bhattacharyya(s1[i], s2[j]) for i in 1:length(s1), j in 1:length(s2)]
end

function dtw_cost_matrix2(seq1::Array{kmerDist}, seq2::Array{kmerDist},
    transportcost=1,
    filterkernel = nothing) where T
    # Build the cost matrix
    m = length(seq2)
    n = length(seq1)

    # Initialize first column and first row
    D = pairwise3(seq1,seq2)
    @assert size(D) == (m,n)

    for r=2:m
        D[r,1] += D[r-1,1]
    end
    for c=2:n
        D[1,c] += D[1,c-1]
    end

    # Complete the cost matrix
    for c = 2:n
        for r = 2:m
            best_neighbor_cost = min(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1])
            D[r, c] += best_neighbor_cost
        end
    end

    if filterkernel !== nothing
        D = imfilter(D, filterkernel)
    end

    return D
end

function trackback2(D::AbstractMatrix{T}) where {T<:Number}

    # initialize trackback throught rows/columns
    r, c       = size(D)
    rows, cols = Int[r], Int[c]

    # estimate that we'll need N⋅logN elements
    N  = max(r, c)
    sz = 2 * N
    sizehint!(rows, sz)
    sizehint!(cols, sz)

    # do trackback
    @inbounds while r > 1 && c > 1
        tb, r, c = indmin3(D[r-1, c-1], D[r-1, c], D[r, c-1], r, c)
        push!(rows, r)
        push!(cols, c)
    end
    # Possibly either r>1 or c>1 at this point (but not both).
    # Add the unfinished part of the track to reach [1,1]
    for r = r-1:-1:1
        push!(rows, r)
        push!(cols, 1)
    end
    for c = c-1:-1:1
        push!(rows, 1)
        push!(cols, c)
    end
    return D[end, end], reverse(cols), reverse(rows)
end

@inline function indmin3(a,b,c,i,j)
    if a <= b
        if a <= c
            return 1,i-1,j-1
        else
            return 3,i,j-1
        end
    else
        if b <= c
            return 2,i-1,j
        else
            return 3,i,j-1
        end
    end
end

function detect_change_points(x::Vector{Float32}, kmers::Dict{Any})
    allKeys = keys(kmers)
    stDv = 0
    for key in allKeys
        stDv += parse(Float64, kmers[key][2])
    end
    stDv /= length(kmers) # Average standard deviation of the kmers list, also defines the threshold
    changePoints = [] # list of changepoints detected
    
    for i in 1:length(x)-1 
        if (abs(x[i] - x[i+1]) > stDv)
            push!(changePoints, i)
        end
    end

    return changePoints
end

