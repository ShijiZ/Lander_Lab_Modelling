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
#print(stats.kstest(p503X_radius, 'norm'))
print(stats.kstest(p503X_radius, 'uniform'))

#plt.figure(1)
#plt.hist(p503X_radius, bins = 100, linewidth = 1,
#         edgecolor = 'black', color = 'cornflowerblue')
#plt.xlabel('Radius')
#plt.ylabel('Frequency')
#plt.title('Nevus radii Distribution of P50 mice treated with 3X (75mg/mL) tamoxifen')

#fitting the data
plt.hist(p503X_radius, bins = 50, linewidth = 1, density=True,
         edgecolor = 'black', color = 'cornflowerblue')
plt.xlabel('Radius ($\mu$m)')
plt.ylabel('Frequency')
plt.title('Nevus radii Distribution of P50 mice treated with 3X (75mg/mL) tamoxifen')

xmin, xmax = plt.xlim()
lnspc = np.linspace(xmin, xmax, len(p503X_radius))

m1, s1, t1 = stats.lognorm.fit(p503X_radius)
pdf_log = stats.lognorm.pdf(lnspc, m1, s1, t1)
plt.plot(lnspc, pdf_log, label = "Lognorm fit", color = "g", linewidth=3)
plt.legend()
plt.show()
