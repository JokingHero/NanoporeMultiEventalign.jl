
"
computes the dtw cost and allignment for sets of fast5 data
"
function nanopore_dtw(x::Vector{Float32}, y::Vector{Float32})
    kmers = loadkmers()
    changepoints1 = detect_change_points(x, kmers)
    changepoints2 = detect_change_points(y, kmers)
    kmerdists1 = kmerdist_from_changepoints(x, changepoints1)
    kmerdists2 = kmerdist_from_changepoints(y, changepoints2)

    return dtw_bhattacharyya(kmerdists1, kmerdists2)
end

"""
computes the dtw of 2 Array{kmerDist} using the bhattacharyya distance
"""
function dtw_bhattacharyya(seq1::Array{kmerDist}, seq2::Array{kmerDist})
    D = dtw_cost_matrix_bhattacharyya(seq1,seq2)
    return trackback_bhattacharyya(D)
end


function dtw_cost_matrix_bhattacharyya(seq1::Array{kmerDist}, seq2::Array{kmerDist},
    transportcost=1) where T
    # Build the cost matrix
    m = length(seq2)
    n = length(seq1)

    # Initialize first column and first row
    D = pairwise_bhattacharyya(seq1,seq2)
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


function trackback_bhattacharyya(D::AbstractMatrix{T}) where {T<:Number}

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
