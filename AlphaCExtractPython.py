# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 23:12:12 2020

@author: Kimberly Taylor
Description: Extracts alpha carbons from a pdb file and writes them to a file
"""
import numpy as np
import csv
#import regex

infile = "C:/Users/kmich/Documents/George Mason/Vaisman Lab/MS Project/PDB/Benchmark 5 Files/benchmark5/structures/1e96_r_u.pdb"
outfile = "C:/Users/kmich/Documents/George Mason/Vaisman Lab/MS Project/PDB/Calpha Extracts/1e96_r_u.pdb"

AllAlphaC = []
count = 0
with open(infile) as f:
    csv_reader = csv.reader(f, delimiter = "\n")
    for row in f: 
        if (row[0:4] == "ATOM") and (row[13:15] == "CA"):
            AllAlphaC = np.append(AllAlphaC,row)
            count += 1

f.close()

g = open(outfile, 'w')

for i in range(len(AllAlphaC)):
    g.write(AllAlphaC[i])

g.close()
