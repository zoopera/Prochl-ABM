#!/usr/bin/env python

import os

ddir = "./aligned_genes"
for filename in os.listdir(ddir):
    if filename.endswith(".fasta"):
	# get column length
	f = open("%s/%s"%(ddir,filename))
	f.readline()
	col = len(f.readline())-1
	f.close()
	# get number of rows (sequences)
        f = open("%s/%s"%(ddir,filename))
	line = f.readline()
	row = 0
	buff = ""
	while line:
	    buff += line[1:]
	    buff += f.readline()
	    row += 1
	    line = f.readline()
	f.close()
	# write to .SEQ file
	fout = open("%s/%s"%(ddir, filename.split(".")[0]+".SEQ"), "w")
	fout.write("%-5d%d\n%s" % (row, col, buff))
	fout.close()
