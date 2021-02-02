
For example data treated sample is extra complicated as some of the bases are methylated and will not conform to the normal distributions. 
Control sample is normal, reference.fa contains reference sequence of our example datasets. Fastq files are basecalled sequences of corresponding fast5 reads.

Fast5 format can be read with H5 libraries, reads from the nanopore have to be normalized before procesing.
PicoAmp normaliation is as follows:

scaling = channel_info["range"] / channel_info["digitisation"]
offset = int(channel_info["offset"]
normalized_data = np.array(scaling * (raw + offset), dtype = np.float32)

each fast5 has its own channel_info and this has to be adjusted for each file.

Reference.fa contains sequence for the expected kmers, combined with the models/r9.4_70bps.u_to_t_rna.5mer.template.model
you can simulate fake reads using level_mean as current mean and level_stdv as current standard deviation.
