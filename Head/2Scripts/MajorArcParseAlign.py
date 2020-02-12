#!/usr/local/bin/python3

import re

score_file = open("score_out", "wt")

for row_index in range(0, 100):
    for column_index in range(0, 100):
        output_file = "_".join([str(row_index), str(column_index)])
        output_path = "/".join(["output", output_file])
        handle = open(output_path, 'rt').read()
        contents = handle.splitlines()
        contents = " ".join(contents)
        score_start = contents.rfind("Score")
        score_string = contents[score_start : score_start + 10]
        temp = re.findall(r'\d+', score_string)
        score_number = list(map(int, temp))[0]
        score_file.write(str(score_number) + "\n")


score_file.close()
