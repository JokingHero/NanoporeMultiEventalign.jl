# Here put file load/file write functions
function loadfast5(path)

    hfile = h5open(path, "r")
    channel_info = hfile["UniqueGlobalKey/channel_id"]
    read_number = "Raw/Reads/" * keys(hfile["Raw/Reads"])[1]
    raw_signal = hfile[read_number]

    # Store all the values used for normalization
    var_range = read(attributes(channel_info)["range"])
    digitisation = read(attributes(channel_info)["digitisation"])
    offset = Int(read(attributes(channel_info)["offset"]))
    raw = read(raw_signal["Signal"]) #Int16 array

    # Normalization
    normalized_data = []
    scaling = var_range / digitisation
    for i in raw
        x = convert(Float32, (scaling * (i + offset)))
        push!(normalized_data, x)
    end
    return normalized_data

end
