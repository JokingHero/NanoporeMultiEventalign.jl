# Here put file load/file write functions
function loadfast5(path)

    hfile = h5open(path, "r")
    channel_info = hfile["UniqueGlobalKey/channel_id"]
    raw_signal = hfile["Raw/Reads/Read_242"]

    # Store all the values used for normalization
    range = read(attributes(channel_info)["range"])
    digitisation = read(attributes(channel_info)["digitisation"])
    offset = Int(read(attributes(channel_info)["offset"]))
    raw = read(raw_signal["Signal"]) #Int16 array

    # Normalization
    normalized_data = []
    scaling = range / digitisation
    for i in raw
        x = scaling * (i + offset)
        push!(normalized_data, x)
    end
    return normalized_data

end
