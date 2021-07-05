#! /Users/femke/miniconda3/bin/python

"""
This script is meant to calculate the nematic order parameter of a system of
rod-like particles based on particle coordinates in a LAMMPS data-file.
"""

import sys
import pandas as pd
import numpy as np
from numpy import linalg as LA

#Assign input values
fp = sys.argv[1]

print("Opening and converting data file....")
f = open(fp, "r")
atoms = f.readlines()[2]
bad_chars = ["a", "t", "o", "m", "s"]
for i in bad_chars:
    atoms = atoms.replace(i, " ")
atoms = int(atoms)

data = pd.read_csv(fp, delimiter="\s+", header=None, skiprows=27, nrows=atoms)
data.columns = ["id", "mol", "type", "x", "y", "z", "nx", "ny", "nz"]
data.set_index("id", inplace=True)
data = data.sort_values("id", axis=0, ascending=True)

CL = data.value_counts(subset="mol")
CL = np.max(CL)

i = 0
j = len(data.index)
while i < j:
    if i + 1 == data.index[i]:
        i += 1
    else:
        molecule = data.iloc[i, 0]
        data = data[data.mol != molecule]
        next = data.loc[data["mol"] == molecule+1]
        i = next.index[0]
        j -= 15

print("Calculating....")
data = data.reset_index()
i = 0
j = len(data.index)
sum = np.zeros((3, 3))
#Iterate over all molecules
while i < j:
    #Calculate the displacements in every direction
    dx = (data.iloc[i+CL-1, 3]) - (data.iloc[i, 3])
    dy = (data.iloc[i+CL-1, 4]) - (data.iloc[i, 4])
    dz = (data.iloc[i+CL-1, 5]) - (data.iloc[i, 5])
    #Calculate the magnitude of the vector
    mag = np.sqrt((dx**2)+(dy**2)+(dz**2))
    #Divide all vector components by the magnitude
    x = dx/mag
    y = dy/mag
    z = dz/mag
    #Calculate tensor components
    xx = ((3/2)*(x*x))-(1/2)
    yy = ((3/2)*(y*y))-(1/2)
    zz = ((3/2)*(z*z))-(1/2)
    xy = ((3/2)*(x*y))
    xz = ((3/2)*(x*z))
    yz = ((3/2)*(y*z))
    #Put the calculated components in an array and add it to the sum for all molecules
    add = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
    sum = np.add(sum, add)
    #Increase i for the next molecule
    i += CL

#Calculate Q and determine the eigenvalues(p) and the corresponding eigenvectors(d)
NC = j/CL #Number of chains in the system
Q = sum/NC
p, d = LA.eig(Q)
print("The eigenvalues are: \n", p)
print("The corresponding eigenvectors are: \n", d)
