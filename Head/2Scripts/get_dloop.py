#!/usr/bin/env python
from Bio import SeqIO
import sys
import os

"""
Notes:
D-loop same as Control region запись вида 

D-loop          15665..16724
                     /note="control region"
"""

outfile = open("/Users/xtinaushakova/WB/tRNA/mito-trna/dloops_control_regions.txt", "wt")
out_stats = open("/Users/xtinaushakova/WB/tRNA/mito-trna/dloop_stats.txt", "wt")

out_stats_header = "\t".join([
	"Species", 
	"D-loop count", 
	"Control region count", 
	"Notes"])
out_stats.write(out_stats_header + "\n")

outfile_header = "\t".join([
	"Species", 
	"Feature name", 
	"Feature location", 
	"Strand", 
	"Sequence", 
	"Notes"])
outfile.write(outfile_header + "\n")


infile = open("/Users/xtinaushakova/WB/tRNA/mito-trna/Body/1Raw/1Source_Genbank/source.gb")

for rec in SeqIO.parse(infile, 'genbank'):
	species = "_".join(rec.annotations['organism'].split())
	control_count = 0
	dloop_count = 0
	same = 0
	for feature in rec.features:
		if feature.type == "D-loop":
			dloop_sequence = str(feature.extract(rec.seq)).upper()
			dloop_location = str(feature.location)
			dloop_strand = str(feature.strand)
			dloop_count += 1
			outfile.write("\t".join([species, "dloop", dloop_location, dloop_strand, dloop_sequence, "", "\n"]))
			if 'notes' in feature.qualifiers and "control region" in feature.qualifiers['note']:
				control += 1
				control_sequence = dloop_sequence
				control_location = dloop_location
				control_strand = dloop_strand
				control_count += 1
				same += 1
				outfile.write("\t".join([species, "control region", control_location, control_strand, control_sequence, "control same as dloop", "\n"]))
		if feature.type == "misc_feature" and 'note' in feature.qualifiers and "control region" in feature.qualifiers['note']:
				control_sequence = str(feature.extract(rec.seq)).upper()
				control_location = str(feature.location)
				control_strand = str(feature.strand)
				control_count += 1
				outfile.write("\t".join([species, "control region", control_location, control_strand, control_sequence, "control annotated as misc_feature", "\n"]))	
	out_stats.write("\t".join([species, str(dloop_count), str(control_count)]))
	if same > 0:
		out_stats.write("\t" + "control same as dloop" + "\n") 
	else: 
		out_stats.write("\t" + "" + "\n")

infile.close()
outfile.close()
out_stats.close()