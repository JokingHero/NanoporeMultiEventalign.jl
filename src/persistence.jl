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

    # Store all the values used for normalization
    var_range = read(attributes(channel_info)["range"])
    digitisation = read(attributes(channel_info)["digitisation"])
    offset = Int(read(attributes(channel_info)["offset"]))
    raw = read(raw_signal["Signal"]) #Int16 array

    # Normalization
    normalized_data = []
    scaling = var_range / digitisation
    normalized_data = map(i -> convert(Float32, (scaling * (i + offset))), raw)
    return normalized_data

end
