# Here put file load/file write functions
function loadfasta(path::String)
    f = FASTA.Reader(open(path,"r"))
    fastaList = []
    for r in f
          push!(fastaList, FASTA.sequence(r))
    end
    close(f)
    return fastaList
end

"""
Makes kmers dictionary that stores mean and stdvt for each kmer.
Access syntax: kmers['GCTGT'][2] gets stdvt of the kmer
"""
function loadkmers(model_path = "models/r9.4_70bps.u_to_t_rna.5mer.template.model")
    all_lines = readlines(model_path)
    kmers = Dict() # Dictionary to read from all the kmers data
    line_counter = 0 # used to skip first 7 lines of the file
    for line in all_lines
        if (line_counter > 6)
            name, mean, stdv = split(line, "\t")
            kmers[name] = mean, stdv
        end
        line_counter = line_counter + 1
    end
    return kmers
end

"""
Gets the corosponding mean and stdv from a fasta sequense using kmer data.
Takes in fasta data loaded from loadfasta and kmers loaded from loadkmers.
Returns a list of tuples where the values are {mean, stdv}.
"""
function fasta_to_kmer_valus(fastadata, kmers)
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
