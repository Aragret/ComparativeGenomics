#!/usr/bin/env python
from Bio import SeqIO
import sys
import os

"""
	This script extracts all annotated dloops from the files you explicitly hand it as an argument. You can either hand it a single genbank file or a whole folder, containing several such files.
	
	usage: python get_dloop.py [input_file_or_dir] [output_file]
"""

# Working function
def get_dloop(filename):
	"""
		This function accepts a genbank file, parsesit, extracts the dloop and returns it as a string of type >name_coords newline sequence
	"""
	dloops = ""
	for rec in SeqIO.parse(filename, 'genbank'):
		species = "_".join(rec.annotations['organism'].split())
		for feature in rec.features:
			if feature.type == "D-loop":
				dloop = (">" + "\t".join([species, str(feature.location)]) + "\n" + str(feature.extract(rec.seq)) + "\n")
				dloops += dloop
	return dloops

#
if len(sys.argv) > 3:
	sys.exit("Too many arguments! Check the doc string for usage info.")
elif len(sys.argv) < 2:
	sys.exit("No input file or dir provided! Check the doc string for usage info.")

user_input = sys.argv[1]

# Check if user specified file or dir as input
dloops = ""
if os.path.isdir(user_input):
	os.chdir(user_input)
	#dloop list of lists
	for filename in os.listdir(user_input):
		dloops += get_dloop(user_input)
	#	dloop = get_dloop(file)
	#	outfile.write(dloop + "\n")
#		print("Здесь мог бы быть ваш длуп")
	#flatten list
elif os.path.isfile(user_input):
	dloops = get_dloop(user_input)

# Check if output file specified
if len(sys.argv) > 2:
	user_output = sys.argv[2]
	user_output = open(user_output, "wt")
	user_output.write(dloops)
	user_output.close()
else:
	print(dloops)