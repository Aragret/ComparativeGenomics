import os
import re
fout = open("OrlovRepeatsInDLoop.txt", 'w')
print('Species', 'Length', 'Number', 'Standard', 'Symm.', 'Compl.', 'Inverted', sep='\t', end='\n', file=fout)


def create_file(filename, data):
    input_filename = filename + ".input"
    output_filename = input_filename + '.out4out'
    with open(input_filename, 'wt') as infile:
        infile.write(data)
    with open(input_filename, 'rt') as infile:
        os.system('OrlovPopadinRepeats2010.exe -m4 -y0 -x10 -d -s -c -i {0}'.format(input_filename))
        infile.close()
        os.remove(input_filename)

    with open(output_filename, 'r') as outfile:
        for line in outfile:
            if re.search("\d{2,3}\t+[1-9]{1}\d*(\t+\d+){4}", line) != None:
                repeats = re.sub('\t\t', '\t', line)
                print(filename, repeats, sep='\t', end='', file=fout)


# "^\d{1,3}\s+[1-9]\d*(\s\d+){4}$"

with open('dloops.fa', 'r') as file:
    for num, line in enumerate(file):
        if num % 2 == 0:
            title = line[1:-1]
            title = title.replace('control region', 'cr')
        else:
            create_file(title, line)

fout.close()
