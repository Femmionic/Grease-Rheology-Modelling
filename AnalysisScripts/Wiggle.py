#! /Users/femke/miniconda3/bin/python

'''This script can be used to extract data from a LAMMPS datafile and analyze it
to extract a storage and loss modulus from it.'''

import numpy as np
import scipy as sp
import scipy.fftpack
import scipy.optimize as optimize
import pandas as pd
import matplotlib.pyplot as plt

file = open("log.wiggle", "r")
lines = file.readlines()
lines = lines[20162:45163]
file.close()

new_file = open("LogWiggle.csv", "w")
new_file.writelines("Step Temp PotEng E_bond TotEng Press Pxy\n")
new_file.writelines(lines)
new_file.close()

data = pd.read_csv("LogWiggle.csv", delim_whitespace=True, usecols=[0,6])

time = []
for item in data["Step"]:
    t = item * 0.005
    time.append(t)
data["Time"] = time

averages = pd.DataFrame(columns=["Time", "Batch_avg"])
i = -390
j = 0
while i < 249600:
    i += 400
    j = i + 390
    batch = data[data["Step"].isin(np.linspace(i, j, 40))]
    avg = batch["Pxy"].mean()
    add = pd.DataFrame([[j*0.005, avg]], columns=["Time", "Batch_avg"])
    averages = averages.append(add)
averages.to_csv("BatchAverage.csv")

time_drop = 750 #Up to which time to drop
data = data.drop(np.arange(0, time_drop/0.05, 1), axis = 0)
averages = averages[averages.Time >= time_drop]

tdata = data["Time"]
Pdata = data["Pxy"]
t_batch = averages["Time"]
P_batch = averages["Batch_avg"]

def stress(t, A, B, v):
    return A*np.sin(2*np.pi*v*t)+B*np.cos(2*np.pi*v*t)

guess_A = (averages["Batch_avg"].max()-averages["Batch_avg"].min())/2
guess_B = guess_A
guess_v = 1/100
guess = np.array([guess_A, guess_B, guess_v])

par, cov = optimize.curve_fit(stress, t_batch, P_batch, guess, bounds=(0.0, np.inf))
print("Sigma':", par[0], "Sigma'':", par[1], "Frequency:", par[2])
Gp = str(par[0])
Gpp = str(par[1])
Freq = str(par[2])

output = open("Calculated.txt", "w")
output.writelines("Sigma': " + Gp + "\n")
output.writelines("Sigma'': " + Gpp + "\n")
output.writelines("Frequency: " + Freq + "\n")
output.close()

plt.scatter(tdata, Pdata, facecolors="none", edgecolors="black", alpha=0.1, label="Instantaneous")
plt.scatter(t_batch, P_batch, facecolors="black", edgecolors="black", label="Batch")
plt.plot(tdata, stress(tdata, par[0], par[1], par[2]), color="Red", label="Fit")
plt.xlim(time_drop, 1250)
plt.xticks([750, 850, 950, 1050, 1150, 1250])
plt.xlabel("Time (\u03C4)")
plt.ylabel(r"$Stress_{xy} (K_{B}T\sigma^{-3})$")
plt.legend()
plt.savefig("StressFit.png", dpi=300)
