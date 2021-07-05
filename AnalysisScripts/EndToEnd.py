"""
This script is meant to calculate the persistence length of semi-flexible rods
based on a LAMMPS data-file.
"""

import sys
import pandas as pd
import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import datetime

start_time = datetime.datetime.now()

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
        i = int(next.index[0])
        j -= 15

print("Converting coordinates....")
data = data.reset_index()
i = 0
j = len(data.index)
#Iterate over all molecules
while i < j:
    data.iloc[i, 3] = (data.iloc[i, 3])+((data.iloc[i, 6]*(5*CL)))
    data.iloc[i, 4] = (data.iloc[i, 4])+((data.iloc[i, 7]*(5*CL)))
    data.iloc[i, 5] = (data.iloc[i, 5])+((data.iloc[i, 8]*(5*CL)))
    i += 1
data = data.drop(["nx", "ny", "nz"], axis=1)

print("Calculating....")
i = 0
j = len(data.index)
EtE = 0
while i < j:
    x = (data.iloc[i+CL-1, 3]) - (data.iloc[i, 3])
    y = (data.iloc[i+CL-1, 4]) - (data.iloc[i, 4])
    z = (data.iloc[i+CL-1, 5]) - (data.iloc[i, 5])
    Mag = np.sqrt((x**2)+(y**2)+(z**2))
    EtE += Mag
    i += CL
EtE = EtE**2
print("Average end-to-end distance = ", EtE)
R2 = EtE/(j/CL)
L = float(CL-1)

print("Fitting the data....")
def End_to_end(L, P):
    return (2*L*P)*(1-(P/L)*(1-np.exp(-L/P)))
pars, cov = curve_fit(f=End_to_end , xdata=L, ydata=R2, p0 = [1], bounds=(0, np.inf))
print("The calculated persistence length is: ", *pars)

print("Runtime: ", datetime.datetime.now()-start_time)
