#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy import stats
import numpy as np

#import all files for same experiment p50 3X
sample32 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1232_tamoxifen3x_P50tab.txt')
sample33 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1233_tamoxifen3x_P50tab.txt')
sample35 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1235_tamoxifen3x_P50tab.txt')
sample36 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1236_tamoxifen3x_P50tab.txt')
sample37 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1237_tamoxifen3x_P50tab.txt')
sample07 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1307_tamoxifen3x_P50tab.txt')
sample09 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1309_tamoxifen3x_P50tab.txt')
sample10 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1310_tamoxifen3x_P50tab.txt')
sample22 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1322_tamoxifen3x_P50tab.txt')
sample30 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/1330_tamoxifen3x_P50tab.txt')
T296 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/T296_PBS_P50_3x.txt')
T306 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/T306_pbs_P50_3x.txt')
T311 = pd.read_table('/Users/shijizhao/Documents/UCI/lab rotation/Lander Lab/Rolando/T311_pbs_P50_3x.txt')

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

CellCycleLimit = 30
r = 3.5 # Radius of a single melanocyte
N = 10 # Number of steps
ProbRange = np.arange(0.4, 1, 0.01) # Probability of transformation
dValueList = []
PercentCellFinalStateList = []
dValue_Data = pd.DataFrame()
PercentCellFinalState_Data = pd.DataFrame()

for n in range(N):
    dValueList.append([])
    PercentCellFinalStateList.append([])
    print ''
    print 'n is', n
    for p in ProbRange:
        FinalRadius = np.array([])
        CellFinalStateList = np.array([])
        for i in range(5000): #Total number of simulations
            CellFinalState = True # Assume all cells goes to final state
            t = 0
            cells = [1]
            for i in range(n):
                cells.append(0) # Initialize the cell in each state

            while sum(cells)-cells[n] != 0: # There are still cells not in the final state
                t += 1
                cells = moveCells(cells, p, n)
                if t > CellCycleLimit: # The cell cycle reaches the limit
                    break

            if sum(cells)-cells[n] != 0:
                CellFinalState = False # There are still cells not in the final state

            FinalRadius = np.append(FinalRadius, r*(np.sum(cells))**(1./3.))
            CellFinalStateList = np.append(CellFinalStateList, CellFinalState)
            #print(cells)
            #print(np.sum(cells))

        dValue = list(stats.ks_2samp(FinalRadius, p503X_radius))[0]
        dValueList[n].append(dValue)
        if dValue < 0.1:
            print 'p is', p
            print dValue
        PercentCellFinalStateList[n].append(np.mean(CellFinalStateList))

dValue_Data['Probability'] = ProbRange
PercentCellFinalState_Data['Probability'] = ProbRange
for n in range(N):
    dValue_Data['n = %d'%n] = dValueList[n]
    PercentCellFinalState_Data['n = %d'%n] = PercentCellFinalStateList[n]
dValue_Data.to_csv('dValue_Data.csv', sep='\t')
PercentCellFinalState_Data.to_csv('Percentage_Data.csv', sep='\t')

plt.figure(1)
for n in range(N):
    plt.plot(ProbRange, dValueList[n], label='n = %d'%n)
plt.xlabel('Probability of Transformation p')
plt.ylabel('d-Value')
plt.title('The d-value Curve against the Transformation Probability p')
plt.legend()

plt.figure(2)
for n in range(N):
    plt.plot(ProbRange, PercentCellFinalStateList[n], label='n = %d'%n)
plt.xlabel('Probability of Transformation p')
plt.ylabel('Percentage of Cells Goes to Final State')
plt.title('The Percentage Curve against the Transformation Probability p')
plt.legend()

plt.show()
