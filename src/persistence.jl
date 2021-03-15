
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
Normalizes and loads nanopore fast5 file.
Takes string that describes full path to the fast5 file as input
and returns the normalized array of type Float32
"""
function loadnanoporefast5(path::String)

    hfile = h5open(path, "r")
    channel_info = hfile["UniqueGlobalKey/channel_id"]
    read_number = "Raw/Reads/" * keys(hfile["Raw/Reads"])[1]
    raw_signal = hfile[read_number]

    # Parameters used for normalization
    var_range = read(attributes(channel_info)["range"])
    digitisation = read(attributes(channel_info)["digitisation"])
    offset = Int(read(attributes(channel_info)["offset"]))
    raw = read(raw_signal["Signal"]) #Int16 array

    # Normalization
    scaling = var_range / digitisation
    normalized_data = map(i -> convert(Float32, (scaling * (i + offset))), raw)
    return normalized_data
end

"
Reads fasta file from `path` and returns a vector of all of the sequences.
"
function loadfasta(path::String)
    f = FASTA.Reader(open(path,"r"))
    seq = []
    for r in f
          push!(seq, FASTA.sequence(r))
    end
    close(f)
    return seq
end
