
"
Probable input: Vector, Vector, distance
"
function dtw_test(fakevar::Bool)
    print("Hello World!")
    return fakevar
end

"
runs the dtw algorthm on the 2 data sets then gets which values from that data
the dtw alignements points to
"
function dtw_stretched(data1, data2)
    list1 = []
    list2 = []
    a,dtw1,dtw2 = dtw(data1,data2)
    for i in dtw1
        push!(list1,data1[i])
    end
    for i in dtw2
        push!(list2,data2[i])
    end
    return list1, list2
end


