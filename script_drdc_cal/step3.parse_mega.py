#!/usr/bin/env python

# This script parses the MEGA-CC output for Ts/Tv rate of each gene family.
# Suppose the *.meg files are under ./mega_out. The output will be "TsTv_ratio.txt"

import os
import sys

with open("TsTv_ratio.txt", "w") as fout:
  for file in os.listdir("./aligned_genes/M10CC_Out/"):
    if file.endswith(".meg"):
	gene_fam = file.split(".")[0]
	# get Ts/Tv rate
        f = open("./aligned_genes/M10CC_Out/" + file)
        line = f.readline()
	iwatch = 0
        while not (line.startswith("[2]  ")):
	    if iwatch > 100000:
	  	sys.exit("Watchdog exited at %s (no '[2]'can be detected)" % file)
   	    line = f.readline()
	    iwatch += 1
	f.close()
        TsTv = line.split("  ")[1]
        # write to file
        fout.write("%s\t%s\n" % (gene_fam, TsTv))

#os.system("ls ./mega_out/*.meg | head -5| xargs mv --target .")
#os.system("rm ./mega_out/*.meg")
#os.system("mv *.meg mega_out")
