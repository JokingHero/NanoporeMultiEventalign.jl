
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
    minVal::Float32
    maxVal::Float32
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
function multi_bhattacharyya(data::Array{kmerdist},
     n::Union{Nothing, Int} = nothing)
     if n === nothing
         lensum = 0
         for d in data
             lensum += length(d.originaldata)
         end
         n = cld(lensum,length(data))
     end
     #findes the maximum and minimum values of all the data
     mini::Float32 = data[1].minVal
     maxi::Float32 = data[1].maxVal
     for d in 2:length(data)
        @inbounds if (mini > data[d].minVal) mini = data[d].minVal end
        @inbounds if (maxi < data[d].maxVal) maxi = data[d].maxVal end
     end
     value = 0;
     partitiondelta::Float32 = (maxi - mini)/n
     for i in 1:n
        part_start::Float32 = mini + partitiondelta * (i-1)
        part_end::Float32 = mini + partitiondelta * i
        part_count = 1
        for d in data
            #counts all the data that is in the interval part_start-part_end
            part_count *= countbetween(part_start, part_end, maxi, d.originaldata)
        end
        # adds the nth-root of the product
        value += part_count^(1/length(data))
     end
     return value
end

"""
found out the built in count function added a lot more
time then it needed to, so this is simple replacement
"""
function countbetween(part_start::Float32,
    part_end::Float32,maxi::Float32,data::Vector{Float32})
    value = 0
    for d in data
        value += (d >= part_start && d < part_end) || d === maxi
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

function multi_pairwise_bhattacharyya(seq::Array{Array{kmerdist}}, sizes::Array{Int64})
    costmat = zeros(sizes...)
    searchpos = convert(Array{Int64},ones(length(seq)))
    searchdata = Array{kmerdist}(undef, length(sizes))
    while true
        for i in 1:length(sizes)
            searchdata[i] = seq[i][searchpos[i]]
        end
        #computes the cost at that point
        costmat[searchpos...] = multi_bhattacharyya(searchdata)

        if searchpos == sizes break end
        for i in 1:length(sizes)
            searchpos[i]+=1
            if searchpos[i] <= sizes[i] break end
            searchpos[i] = 1
        end
    end
    return costmat
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
        push!(kmerdists, kmerdist(average, stdv, min(datalist...), max(datalist...), datalist))
    end
    return kmerdists
end
