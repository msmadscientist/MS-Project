# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 23:53:26 2019

@author: Kimberly Taylor
Title: SurfaceResiduesAllAtoms
Description: Finds surface residues using all atoms in residues using the 
             QUAD method
             
12/18/2019 - this works.  Best coordination with PyMol seems to be with 12 A depth for surface and
5 A minimum distance for interface.
"""

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#Above line must be included even though it is not technically used

import numpy as np
import csv
import matplotlib.pyplot as mpl
import collections
from scipy import optimize

#Calculates value of the ellipsoid at point x, y, z
def ellipsoid (x,y,z,xc,yc,zc,a,b,c):
    return ((x-xc)/a)**2 + ((y-yc)/b)**2 + ((z-zc)/c)**2 

#Calculate distance between two points
def distance (x1, y1, z1, x2, y2, z2):
    return ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**(1/2)

#Calculate center of ellipse
def CalcCenter(xc, yc, zc):
#Calculate the center of the ellipsoid using estimated radii for three axes
   return (((xcoords - xc))**2/(xr**2) + ((ycoords - yc))**2/(yr**2) + ((zcoords - zc))**2/(zr**2))

def f_2 (c):
#Error function for calculation of the center of the ellipsoid
    Ci = CalcCenter(*c)
    return Ci - Ci.mean()

def CleanArray (a, b):
    for i in a:
        index = np.where(b == i)    #Get the index in b for value i
        b = np.delete(b,index)
            
    return b

pdb = "2J7P"
infile = "C:/Users/kmich/Documents/George Mason/Vaisman Lab/MS Project/PDB/Original PDB/"+pdb+".pdb"
outfile1 = "C:/Users/kmich/Documents/George Mason/Vaisman Lab/MS Project/PDB/QuAD Results/"+pdb+"_allatom_surface5.csv"
outfile2 = "C:/Users/kmich/Documents/George Mason/Vaisman Lab/MS Project/PDB/QuAD Results/"+pdb+"_allatom_interface5.csv"
chains = np.array(['A','D'])    #Selected chains to calculate surface and interface residues


#Read in coordinates for alpha carbons only.  Save them in a dictionary with keys
#for each chain in the pdb file
allresidues = collections.defaultdict(dict) #All residues in protein
surfaceresidues = collections.defaultdict(dict) #All surface residues

linenum = 0

#Define list of amino acids
AllAminos = np.array(['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'])

#Read coordinates from file into dictionary, keyed by subunit
with open(infile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = "\n")
    for row in csv_file:      
# New code using column lengths should solve column issues
        if (row[0:4] == "ATOM"):# and (row[17:20] in AllAminos):
            if len(allresidues[row[21]]) == 0: # Assignment for first instance of key
                allresidues[row[21]] = np.asarray((float(row[22:26]),float(row[30:38]),float(row[38:46]),float(row[46:54])))
            else: #If key exists, add coordinates to key entry
                allresidues[row[21]] = np.vstack((allresidues[row[21]],np.asarray((float(row[22:26]),float(row[30:38]),float(row[38:46]),float(row[46:54])))))

#Old code for splitting on whitespace - doesn't always work
#        row = row.split()
#        if (row[0] == "ATOM") and (row[3] in AllAminos):
#            if len(allresidues[row[4]]) == 0: # Assignment for first instance of key
#                allresidues[row[4]] = np.asarray((float(row[5]),float(row[6]),float(row[7]),float(row[8])))
#            else: #If key exists, add coordinates to key entry
#                allresidues[row[4]] = np.vstack((allresidues[row[4]],np.asarray((float(row[5]),float(row[6]),float(row[7]),float(row[8])))))


csv_file.close()

keycombo = np.asarray([*chains]) #Used to obtain combinatorial of keys
radii = collections.defaultdict(dict)
centers = collections.defaultdict(dict)
depth = 12 #Use 12 angstroms as the average depth of the surface layer of residues

for key in chains:
    resnums = allresidues[key][:,0]
    xcoords = allresidues[key][:,1]
    ycoords = allresidues[key][:,2]
    zcoords = allresidues[key][:,3]
    
#Get estimated radii for each axis 
    xr = abs(np.max(xcoords) - np.min(xcoords))/2
    yr = abs(np.max(ycoords) - np.min(ycoords))/2
    zr = abs(np.max(zcoords) - np.min(zcoords))/2
    radii[key] = np.asarray((xr,yr,zr))

#Perform fitting for each chain in the pdb file    
    EstCenter = (15,15,15)
    centers[key], ier = optimize.leastsq(f_2,(EstCenter))

#Calculate metrics and report results
    xc_2, yc_2, zc_2 = centers[key]
    Ri_2       = CalcCenter(xc_2,yc_2,zc_2)
    R_2        = Ri_2.mean()
    residu_2   = sum((Ri_2 - R_2)**2)

#Find surface residues    
    for n in range(len(xcoords)):
        if (ellipsoid(xcoords[n],ycoords[n],zcoords[n],xc_2,yc_2,zc_2,xr,yr,zr) <= 1) and\
        ((ellipsoid(xcoords[n],ycoords[n],zcoords[n],xc_2,yc_2,zc_2,xr-depth,yr-depth,zr-depth) >= 1)):
            if key not in surfaceresidues:
                surfaceresidues[key] = np.asarray(resnums[n])
            else:
                surfaceresidues[key] = np.append(surfaceresidues[key],np.asarray(resnums[n]))

    surfaceresidues[key] = np.unique(surfaceresidues[key])
    


#Use surface residues to find interface residues
mindist = 5 # Minimum recommended distance is 5 Angstroms
interresidues = collections.defaultdict(dict) #All surface residues
count = 0
if len(keycombo) > 1:
    for i in range(len(keycombo)-1):    #i is one monomer, j is the other
        for j in range(i+1,len(keycombo)):
            x1 = allresidues[keycombo[i]][:,1]
            y1 = allresidues[keycombo[i]][:,2]
            z1 = allresidues[keycombo[i]][:,3]
            resnum1 = allresidues[keycombo[i]][:,0]
            x2 = allresidues[keycombo[j]][:,1]
            y2 = allresidues[keycombo[j]][:,2]
            z2 = allresidues[keycombo[j]][:,3]
            resnum2 = allresidues[keycombo[j]][:,0]
            for m in range(len(x1)):
                for n in range(len(x2)):
                    if (resnum1[m] in surfaceresidues[keycombo[i]]) and (resnum2[n] in surfaceresidues[keycombo[j]]):
                        count += 1
                        if distance(x1[m],y1[m],z1[m],x2[n],y2[n],z2[n]) < mindist:
                            if keycombo[i] not in interresidues:
                                interresidues[keycombo[i]] = np.asarray(resnum1[m])
                            else:
                                interresidues[keycombo[i]] = np.append(interresidues[keycombo[i]],resnum1[m])
                            if keycombo[j] not in interresidues:
                                interresidues[keycombo[j]] = np.asarray(resnum2[n])
                            else:
                                interresidues[keycombo[j]] = np.append(interresidues[keycombo[j]],resnum2[n])
                            
    for key in interresidues:
        interresidues[key] = np.unique(interresidues[key])

#Remove interface residues from surface residues
        surfaceresidues[key] = CleanArray(interresidues[key],surfaceresidues[key]) 
 
    #Output interface residues
    with open(outfile2, 'w', newline = '') as csv_file2:
        rowwriter = csv.writer(csv_file2, delimiter = ',' )
        for key in interresidues:
            for n in range(len(interresidues[key])):
                rowwriter.writerow((key,interresidues[key][n]))
    
    csv_file2.close()
    
          
                
#Output surface residues
with open(outfile1, 'w', newline = '') as csv_file1:
    rowwriter = csv.writer(csv_file1, delimiter = ',' )
    for key in surfaceresidues:
        for n in range(len(surfaceresidues[key])):
            rowwriter.writerow((key,surfaceresidues[key][n]))

csv_file1.close()    
