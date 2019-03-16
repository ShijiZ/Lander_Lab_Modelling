#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy import stats
import numpy as np

path = '/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Lander_LAB_Modelling/Rolando/'

#import all files for same experiment p50 3X
sample32 = pd.read_table('%s1232_tamoxifen3x_P50tab.txt'%path)
sample33 = pd.read_table('%s1233_tamoxifen3x_P50tab.txt'%path)
sample35 = pd.read_table('%s1235_tamoxifen3x_P50tab.txt'%path)
sample36 = pd.read_table('%s1236_tamoxifen3x_P50tab.txt'%path)
sample37 = pd.read_table('%s1237_tamoxifen3x_P50tab.txt'%path)
sample07 = pd.read_table('%s1307_tamoxifen3x_P50tab.txt'%path)
sample09 = pd.read_table('%s1309_tamoxifen3x_P50tab.txt'%path)
sample10 = pd.read_table('%s1310_tamoxifen3x_P50tab.txt'%path)
sample22 = pd.read_table('%s1322_tamoxifen3x_P50tab.txt'%path)
sample30 = pd.read_table('%s1330_tamoxifen3x_P50tab.txt'%path)
T296 = pd.read_table('%sT296_PBS_P50_3x.txt'%path)
T306 = pd.read_table('%sT306_pbs_P50_3x.txt'%path)
T311 = pd.read_table('%sT311_pbs_P50_3x.txt'%path)

#index the column of interest which in this case was nevus area
sample32 = sample32.iloc[:,1]
sample33 = sample33.iloc[:,1]
sample35 = sample35.iloc[:,1]
sample36 = sample36.iloc[:,1]
sample37 = sample37.iloc[:,1]
sample07 = sample07.iloc[:,1]
sample09 = sample09.iloc[:,1]
sample10 = sample10.iloc[:,1]
sample22 = sample22.iloc[:,1]
sample30 = sample30.iloc[:,1]
T296 = T296.iloc[:,1]
T306 = T306.iloc[:,1]
T311 = T311.iloc[:,1]

#append all files into one dataframe (df)
p503X = sample32.append([sample33, sample35, sample36, sample37, sample07,
                         sample09, sample10, sample22, sample30, T296,
                         T306, T311])

#remove NaN
p503X = p503X.dropna()
p503X = p503X[p503X > 0]

#get the radus
p503X_radius = p503X.apply(lambda x: math.sqrt(x/math.pi))

def moveCells(cells, p, n):

    m = []
    for i in range(n):
        if cells[i] == 0:
            m.append(0)
        else:
            m.append(min(cells[i], np.random.poisson(lam=(p*cells[i]))))

    Newcells = []
    Newcells.append(2*(cells[0]-m[0]))
    for i in range(n-1):
        Newcells.append(2*(cells[i+1]+m[i]-m[i+1]))
    Newcells.append(cells[n]+m[n-1])

    return Newcells

n = 5 # Number of steps
p = 0.74 # Probability of transformation
CellCycleLimit = 30
r = 3.5 # The radius of a single melanocyte

FinalRadius = np.array([])
for i in range(5000): # Total number of simulations
    t = 0
    cells = [1]
    for i in range(n):
        cells.append(0)

    while sum(cells)-cells[n] != 0: # There are still cells not in the final state
        t += 1
        cells = moveCells(cells, p, n)
        if t > CellCycleLimit: # The cell cycle reaches the limit
            break

    FinalRadius = np.append(FinalRadius, r*(np.sum(cells))**(1./3.))
    print(cells)
    print(np.sum(cells))

values1, base1 = np.histogram(p503X_radius, bins=50, density=True)
values2, base2 = np.histogram(FinalRadius, bins=50, density=True)
#evaluate the cumulative
cumulative1 = np.cumsum(values1)
cumulative2 = np.cumsum(values2)
# plot the cumulative function
plt.plot(base1[:-1], cumulative1, c='blue', label='Experimental data')
plt.plot(base2[:-1], cumulative2, c='red', label='Simulated data')
plt.xlabel('Radius ($\mu$m)')
plt.ylabel('Fraction')
plt.title('Comparison of Cumulative Distributions of \n'
          'Experimental and Simulated Data')
plt.legend()

plt.show()
