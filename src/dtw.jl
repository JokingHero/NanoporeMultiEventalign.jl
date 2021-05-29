
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


function multi_nanopore_dtw(x::Vector{Vector{Float32}},
     kmerpath::String = "models/r9.4_70bps.u_to_t_rna.5mer.template.model")
     kmers = loadkmers(kmerpath)
     kmerdists::Array{Array{kmerdist}} = []
     for i in x
        push!(kmerdists, kmerdist_from_changepoints(i, detect_change_points(i, kmers)))
     end
     seqlengths::Array{Int64} = length.(kmerdists)
     D = multi_pairwise_bhattacharyya(kmerdists, seqlengths)
     return multi_trackback_bhattacharyya(D, seqlengths)
end

"""
temporary function to graphically test the multi nanopore dtw function
this only works on 2 sequences 
"""
function multi_nanopore_dtw_plot(x::Vector{Vector{Float32}},
    seperation::Int64 = 20, stepsize::Int64 = 20,
     kmerpath::String = "models/r9.4_70bps.u_to_t_rna.5mer.template.model")
     kmers = loadkmers(kmerpath)
     kmerdists::Array{Array{kmerdist}} = []
     for i in x
        push!(kmerdists, kmerdist_from_changepoints(i, detect_change_points(i, kmers)))
     end
     seqlengths::Array{Int64} = []
     for s in kmerdists
         push!(seqlengths, length(s))
     end
     D = multi_pairwise_bhattacharyya(kmerdists, seqlengths)
     a,b = multi_trackback_bhattacharyya(D, seqlengths)
     o = []
     for i in kmerdists[1]
         push!(o,i.mean+seperation)
     end
     pl = plot(o,thickness_scaling = 2, size = (1500,1000)) #plots mean points for the first graph
     p = []
     for i in kmerdists[2]
         push!(p,i.mean-seperation)
     end
     plot!(p,thickness_scaling = 2) #plots mean points for the second graph

     for i in 1:stepsize:length(b[1])# plots the connecting lines between the grafs
         plot!([b[1][i], b[2][i]], [o[b[1][i]], p[b[2][i]]],color = :black,thickness_scaling = 1, alpha = 0.5, label = false)
     end
     return pl
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

function multi_trackback_bhattacharyya(D::Array{Float64},
     sizes::Vector{Int64})

    # estimate that we'll need N⋅logN elements
    N  = max(sizes...)
    sz = 2 * N
    data = [Int[] for i=1:length(sizes)]
    for i in 1:length(sizes)
        push!(data[i], sizes[i])
    end

    #starting point
    point = copy(sizes)
    nullp = convert(Vector{Int64},ones(length(sizes)))
    # do trackback
    @inbounds while point != nullp
        point = indminn(D,point)
        for i in 1:length(point)
            push!(data[i], point[i])
        end
    end
    for i in 1:length(data)
        reverse!(data[i])
    end
    return D[sizes...], data
end


"""
a function to find the indexes of the minimum value out of
costmat[x-1,y-1,...], costmat[x-1,y,...], costmat[x,y-1,...], ect
at a pos when costmat is an n-dimensional matrix
"""
function indminn(costmat, pos::Array{Int64})
    minval = Inf
    len = length(pos)
    bestpos = Vector{Int64}(undef,len)

    for i in 1:((1 << len)-1)
        newpos = copy(pos)
        for j in 1:len
            newpos[j] -= (0!= i&(1<<(j-1)))
        end
        if !checkbounds(Bool, costmat, newpos... ) continue end
        val = costmat[newpos...]
        if (val <= minval)
            minval = val
            bestpos = newpos
        end
    end
    return bestpos
end
