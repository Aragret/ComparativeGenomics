#!/usr/local/bin/python3

import os
from Bio import SeqIO

# Get sequence file, read as long string
mitoref_file = "1_mitoref_majorarc.fasta"
with open(mitoref_file, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        mitoref_seq = str(record.seq)

for row_index, i in enumerate(range(0, len(mitoref_seq), 100)):
    seq_short = mitoref_seq[i:i + 100]
    filename = "input_" + str(row_index)
    filepath = "/".join(["input", filename])
    input = open(filepath, "wt")
    input.write(seq_short)
    input.close()
    #print(i, filename)
    #print(seq_short)
    #print(len(mitoref_seq))
    for column_index, j in enumerate(range(0, len(mitoref_seq), 100)):
        input_i = "input_" + str(row_index)
        filepath_i = "/".join(["input", input_i])
        input_j = "input_" + str(column_index)
        filepath_j = "/".join(["input", input_j])
        output_file = "_".join([str(row_index), str(column_index)])
        output_path = "/".join(["output", output_file])
        command = "needle -gapopen 10.0 -gapextend 0.5 -asequence " + filepath_i + " -bsequence " + filepath_j + " -outfile " + output_path
        os.system(command)
        # To keep track
        print("row = %d, j = %d" % (row_index, i))
        print("column = %s, n = %s" % (column_index, j))
