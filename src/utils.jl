

function plotdata(path::String)
    kmers = loadkmers()
    data = loadfasta(path)
    dataString = convert(String, data[1])
    plotList = []
    for i in 1:(length(dataString)-4)
        kmerKey = dataString[i:(i+4)]
        #collects all the kmer values of the data into a lits
        push!(plotList,parse(Float64,kmers[kmerKey][1]))
    end
    plot(1:(length(plotList)),plotList, )
    png("plot")#saves the plot as a png named plot
end
