# Here put file load/file write functions
function loadReference(path)
    f = FASTA.Reader(open(path,"r"))
    for r in f
         return FASTA.sequence(r)
    end
end
