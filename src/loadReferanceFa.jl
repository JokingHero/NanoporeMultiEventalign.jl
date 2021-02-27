#a function to load in referance data
function loadReference(reference)
    f = open(string("test/sample_data/",  reference))
    foo = split(read(f,String),'\n')[2]
    p = chop(foo)
    seq = ReferenceSequence(p)
    return seq
end

print(loadReference("reference.fa"))
