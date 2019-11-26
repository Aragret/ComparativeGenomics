import sys
from Bio import AlignIO
import subprocess

print('first argument is alignment, second is refseq fasta')

alignment = AlignIO.read(sys.argv[1], format='clustal')
alignment.sort()

AlignIO.write(alignment, sys.argv[1] + '.sorted', 'clustal')

subprocess.call("seqkit sort --quiet -i {0} > {1}".format(sys.argv[2], sys.argv[2] + '.sorted'), shell=True)

subprocess.call('perl pal2nal.pl -codontable 5 {0} {1} > {2}'.format(sys.argv[1] + '.sorted', sys.argv[2] + '.sorted', sys.argv[2] + '.codons'), shell=True)
