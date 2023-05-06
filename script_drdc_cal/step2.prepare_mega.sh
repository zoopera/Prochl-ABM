#!/bin/bash

# generate list file
find ./aligned_genes/ -name "*fasta" > filelist.txt

# run MEGA-cc
~/software/megacc -a ts2tv.mao -d filelist.txt
