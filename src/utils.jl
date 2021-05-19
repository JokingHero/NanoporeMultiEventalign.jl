
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

struct kmerdist
    mean::Float32
    sd::Float32
    originaldata::Array{Float32}
end


"""
computes the Bhattacharyya distance of two kmerDist x and y
"""
function bhattacharyya(x::kmerdist,y::kmerdist)
    xsd = x.sd
    ysd = y.sd
    if xsd == 0
        xsd = 0.000000001
    end
    if ysd == 0
        ysd = 0.000000001
    end
    arg1 = 1/4*((xsd*xsd)/(ysd*ysd) + (ysd*ysd)/(xsd*xsd) + 2)
    arg2 = 1/4*((x.mean - y.mean)*(x.mean - y.mean)/(xsd*xsd + ysd*ysd))
    return 1/4*log(arg1 + arg2)
end

"""
a function that computes the Bhattacharyya coefficient of
n-distribution
"""
function multi_bhattacharyya(data::Array{Array{Float32}})
    value = 0
    for i in 1:length(data[1])
        d = 1
        for j in 1:length(data)
            d *= data[j][i]
        end
        #takes the n-th root of the product and adds it to the total
        value += d^(1/length(data))
    end
    return value
end

"""
computetes the bhattacharyya distance for every pair in 2 kmerDist arrays
and fills out an array with the result
"""
function pairwise_bhattacharyya(s1::Array{kmerdist},s2::Array{kmerdist})
    return [bhattacharyya(s1[j], s2[i]) for i in 1:length(s2), j in 1:length(s1)]
end


function multi_pairwise_bhattacharyya(seq::Array{Array{Array{Float32}}}, sizes::Array{Int64})
    costmat = zeros(sizes...)
    searchpos = convert(Array{Int64},ones(length(seq)))
    multi_fill_costmatrix(seq, costmat, searchpos)
    return costmat
end

"""
a function that recusivly fills out a cost matrix (costmat)
"""
function multi_fill_costmatrix(seq::Array{Array{Array{Float32}}}, costmat, searchpos)
    #checks if thos point has not already been computed
    if costmat[searchpos...] != 0
        return
    end
    searchdata::Array{Array{Float32}} = []
    for i in 1:length(seq)
        push!(searchdata, seq[i][searchpos[i]])
        #creates a new point to compute
        newsearch = copy(searchpos)
        newsearch[i]+=1
        #checks if this point is inside of the cost matrix
        if checkbounds(Bool,costmat, newsearch...)
            multi_fill_costmatrix(seq, costmat, newsearch)
        end
    end
    #computes the cost at that point
    costmat[searchpos...] = multi_bhattacharyya(searchdata)
end

@inline function indmin3(a,b,c,i,j)
    if a <= b
        if a <= c
            return 1, i-1, j-1
        else
            return 3, i, j-1
        end
    else
        if b <= c
            return 2, i-1, j
        else
            return 3, i, j-1
        end
    end
end

"""
detects and returns a list of change points in fast5 data
using the average kmers standard deviation as a threshold
"""
function detect_change_points(x::Vector{Float32}, kmers::Dict{Any})
    allKeys = keys(kmers)
    stDv = 0
    for key in allKeys
        stDv += parse(Float64, kmers[key][2])
    end
    stDv /= length(kmers) # Average standard deviation of the kmers list, also defines the threshold
    changePoints = [] # list of changepoints detected
    push!(changePoints, 1)
    for i in 1:length(x)-1
        if (abs(x[i] - x[i+1]) > stDv)
            push!(changePoints, i+1)
        end
    end
    push!(changePoints, length(x)+1)
    return changePoints
end



"""
computes the mean and standard deviation for each of the changepoints
in fast5 data
"""
function kmerdist_from_changepoints(x::Vector{Float32}, changepoints::Vector{Any})::Array{kmerdist}
    kmerdists = []
    for i in 2:length(changepoints)
        #computes the average of a changepoint
        average = 0
        datalist = []
        for j in changepoints[i-1]:changepoints[i]-1
            average += x[j]
            push!(datalist,x[j])
        end
        average /= changepoints[i] - changepoints[i-1]
        #computes the standard deviation of a changeoint
        stdv = 0
        for j in changepoints[i-1]:changepoints[i]-1
            stdv += (x[j] - average)*(x[j] - average)
        end
        stdv /= changepoints[i] - changepoints[i-1]
        stdv = sqrt(stdv)
        push!(kmerdists, kmerdist(average, stdv, datalist))
    end
    return kmerdists
end
