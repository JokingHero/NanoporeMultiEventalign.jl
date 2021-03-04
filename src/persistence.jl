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
