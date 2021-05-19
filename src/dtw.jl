
"
computes the dtw cost and allignment for sets of fast5 data
"
function nanopore_dtw(x::Vector{Float32}, y::Vector{Float32},
     kmerpath::String = "models/r9.4_70bps.u_to_t_rna.5mer.template.model")
    kmers = loadkmers(kmerpath)
    changepoints1 = detect_change_points(x, kmers)
    changepoints2 = detect_change_points(y, kmers)
    kmerdists1 = kmerdist_from_changepoints(x, changepoints1)
    kmerdists2 = kmerdist_from_changepoints(y, changepoints2)
    return dtw_bhattacharyya(kmerdists1, kmerdists2)
end


function nanopore_multi_dtw(x::Vector{Vector{Float32}},
     kmerpath::String = "models/r9.4_70bps.u_to_t_rna.5mer.template.model")

    return dtw_multiple_cost_matrix_bhattacharyya(x)
end

"
computes the dtw cost and allignment for sets of fast5 data and plots the result
"
function nanopore_dtw_plot(x::Vector{Float32}, y::Vector{Float32},
     seperation::Int64 = 20, stepsize::Int64 = 20,
     kmerpath::String = "models/r9.4_70bps.u_to_t_rna.5mer.template.model")
    kmers = loadkmers()
    changepoints1 = detect_change_points(x, kmers)
    changepoints2 = detect_change_points(y, kmers)
    kmerdists1 = kmerdist_from_changepoints(x, changepoints1)
    kmerdists2 = kmerdist_from_changepoints(y, changepoints2)
    o = []
    for i in kmerdists1
        push!(o,i.mean+seperation)
    end
    pl = plot(o,thickness_scaling = 2, size = (1500,1000)) #plots mean points for the first graph
    p = []
    for i in kmerdists2
        push!(p,i.mean-seperation)
    end
    plot!(p,thickness_scaling = 2) #plots mean points for the second graph
    a,b,c = dtw_bhattacharyya(kmerdists1, kmerdists2)

    for i in 1:stepsize:length(b) # plots the connecting lines between the grafs
        plot!([b[i], c[i]], [o[b[i]], p[c[i]]],color = :black,thickness_scaling = 1, alpha = 0.5, label = false)
    end
    return pl
end

"""
computes the dtw of 2 Array{kmerDist} using the bhattacharyya distance
"""
function dtw_bhattacharyya(seq1::Array{kmerdist}, seq2::Array{kmerdist})
    D = dtw_cost_matrix_bhattacharyya(seq1,seq2)
    return trackback_bhattacharyya(D)
end

function dtw_multiple_cost_matrix_bhattacharyya(seq::Array{Array{Float32}},
    segmentsize = 100, transportcost=1) where T

    #split the data into segments of size segmentsize
    #note: if the sequences doesn't cleanly split into
    #segments, it will just discard the remaining data
    splitseq::Array{Array{Array{Float32}}} = []
    for s in seq
        i = 1
        segments = []
        while i+segmentsize-1 < length(s)
            push!(segments,s[i:i+segmentsize-1])
            i += segmentsize
        end
        push!(splitseq, segments)
    end

    #stores the length of all the sequenses
    seqlengths::Array{Int64} = []
    for s in splitseq
        if length(s) == 0
            return
        end
        push!(seqlengths, length(s))
    end

    # Build the cost matrix
    D = multi_pairwise_bhattacharyya(splitseq, seqlengths)

    # Initialize first column and first row
    # for i in 1:length(D)
    #     for r=2:size(D[i])[1]
    #         D[i][r, 1] += D[i][r-1, 1]
    #     end
    #     for c=2:size(D[i])[2]
    #         D[i][1, c] += D[i][1, c-1]
    #     end
    # end
    #
    #
    # # Complete the cost matrix
    # for i in 1:length(D)
    #     for c = 2:size(D[i])[2]
    #         for r = 2:size(D[i])[1]
    #             best_neighbor_cost = min(transportcost*D[i][r-1, c], D[i][r-1, c-1], transportcost*D[i][r, c-1])
    #             D[i][r, c] += best_neighbor_cost
    #         end
    #     end
    # end

    return D
end


function dtw_cost_matrix_bhattacharyya(seq1::Array{kmerdist}, seq2::Array{kmerdist},
    transportcost=1) where T
    # Build the cost matrix
    m = length(seq2)
    n = length(seq1)

    # Initialize first column and first row
    D = pairwise_bhattacharyya(seq1,seq2)
    @assert size(D) == (m, n)

    for r=2:m
        D[r, 1] += D[r-1, 1]
    end
    for c=2:n
        D[1, c] += D[1, c-1]
    end

    # Complete the cost matrix
    for c = 2:n
        for r = 2:m
            best_neighbor_cost = min(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1])
            D[r, c] += best_neighbor_cost
        end
    end

    return D
end


function trackback_bhattacharyya(D::AbstractMatrix{T}) where {T<:Number}

    # initialize trackback throught rows/columns
    r, c = size(D)
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
