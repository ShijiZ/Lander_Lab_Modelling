#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 13:40:01 2018

@author: Ruiz
"""

import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy import stats
import numpy as np

#import all files for same experiment p50 3X
sample32 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1232_tamoxifen3x_P50tab.txt')
sample33 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1233_tamoxifen3x_P50tab.txt')
sample35 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1235_tamoxifen3x_P50tab.txt')
sample36 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1236_tamoxifen3x_P50tab.txt')
sample37 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1237_tamoxifen3x_P50tab.txt')
sample07 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1307_tamoxifen3x_P50tab.txt')
sample09 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1309_tamoxifen3x_P50tab.txt')
sample10 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1310_tamoxifen3x_P50tab.txt')
sample22 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1322_tamoxifen3x_P50tab.txt')
sample30 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/1330_tamoxifen3x_P50tab.txt')
T296 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P50_tabs/T296_PBS_P50_3x.txt')
T306 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P50_tabs/T306_pbs_P50_3x.txt')
T311 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P50_tabs/T311_pbs_P50_3x.txt')

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

#make a histogram
# =============================================================================
# plt.hist(p503X_radius, bins = 25, linewidth = 1, edgecolor = 'black',
#          color = 'cornflowerblue')
# plt.xlabel('Nevus area (sq µm)')
# plt.ylabel('Frequency')
# plt.title('Nevus area Distribution of P50 mice treated with 3X (75mg/mL) tamoxifen')
# plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/04162018_nevusAreaDistribution.pdf')
# 
# =============================================================================
plt.hist(p503X_radius, bins = 30, range = (0,65), linewidth = 1,
         edgecolor = 'black', color = 'cornflowerblue')
plt.ylim(0, 55)
plt.xlabel('Radius')
plt.ylabel('Frequency')
plt.title('Nevus radii Distribution of P50 mice treated with 3X (75mg/mL) tamoxifen')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P50_3x/04162018nevusRadiDistribution.pdf')

#fitting ht edata
plt.hist(p503X_radius, bins = 30, linewidth = 1,
         edgecolor = 'black', color = 'cornflowerblue', normed = True)
xt = plt.xticks()[0]
xmin, xmax = min(xt), max(xt)

lnspc = np.linspace(xmin, xmax, len(p503X_radius))

m, s = stats.norm.fit(p503X_radius) 
pdf_g = stats.norm.pdf(lnspc, m, s)  
plt.plot(lnspc, pdf_g, label="Norm", color = "plum", linewidth=3) 

m1, s1, t1 = stats.lognorm.fit(p503X_radius)
pdf_log = stats.lognorm.pdf(lnspc, m1, s1, t1)
plt.plot(lnspc, pdf_log, label = "lognorm", color = "lightgreen", linewidth=3)

plt.legend()

#import all files for same experiment p21 1X
sample39 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_1x/1339_tamoxifen1x_P21tab.txt')
sample41 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_1x/1341_tamoxifen1x_P21tab.txt')
sample44 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_1x/1344_tamoxifen1x_P21tab.txt')
sample50 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_1x/1350_tamoxifen1x_ P21tab.txt')

#index the column of interest which in this case was nevus area
sample39 = sample39.iloc[:,1]
sample41 = sample41.iloc[:,1]
sample44 = sample44.iloc[:,1]
sample50 = sample50.iloc[:,1]

#append all files into one dataframe (df)
p211X = sample39.append([sample41, sample44, sample50])

#remove NaN
p211X = p211X.dropna()
p211X = p211X[p211X > 0]

#make a histogram
plt.hist(df, bins = 25, linewidth = 1, edgecolor = 'black')
plt.xlabel('Nevus area (sq µm)')
plt.ylabel('Frequency')
plt.title('Nevus radii Distribution of P21 mice treated with 1X (25mg/mL) tamoxifen')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P21_1x/nevusAreaDistribution.pdf')

plt.hist(p211X_radius, bins = 30, range = (0,60), linewidth = 1, edgecolor = 'black',
         color = 'cornflowerblue')
plt.ylim(0, 40)
plt.xlabel('Radius (sq µm)')
plt.ylabel('Frequency')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P21_1x/nevusRadiDistribution.pdf')


#import all files for same experiment p21 3X
sample32 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1232_tamoxifen3x_P21tab.txt')
sample33 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1233_tamoxifen3x_P21tab.txt')
sample35 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1235_tamoxifen3x_P21tab.txt')
sample36 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1236_tamoxifen3x_P21tab.txt')
sample37 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1237_tamoxifen3x_P21tab.txt')
sample07 = pd.read_csv('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1307_tamoxifen3x_P21tab.csv')
sample09 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1309_tamoxifen3x_P21tab.txt')
sample10 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1310_tamoxifen3x_P21tab.txt')
sample22 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1322_tamoxifen3x_P21tab.txt')
sample30 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/1330_tamoxifen3x_P21tab.txt')
sample178 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/ATR/T178_WT3x_PBScont_MFtab.txt')
T296 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P21_tabs/T296_PBS_P21_3x.txt')
T305 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P21_tabs/T305_PBS_P21_3x.txt')
T306 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P21_tabs/T306_PBS_P21_3x.txt')
T311 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P21_tabs/T311_PBS_P21_3x.txt')

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
p21 = sample178.iloc[:,1]
T296 = T296.iloc[:,1]
T305 = T305.iloc[:,1]
T306 = T306.iloc[:,1]
T311 = T311.iloc[:,1]

#append all files into one dataframe (df)
p213X = sample32.append([sample33, sample35, sample36, sample37, sample07, 
                         sample09, sample10, sample22, sample30, p21, T296,
                         T305, T306, T311])

#remove NaN
p213X = p213X.dropna()
p213X = p213X[p213X > 0]

#get the radius
p213X_radius = p213X.apply(lambda x: math.sqrt(x/math.pi))

#make a histogram
# =============================================================================
# plt.hist(df, bins = 25, linewidth = 1, edgecolor = 'black')
# plt.xlabel('Nevus area (sq µm)')
# plt.ylabel('Frequency')
# plt.title('Nevus area Distribution of P21 mice treated with 3X (75mg/mL) tamoxifen')
# plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/nevusAreaDistribution.pdf')
# 
# =============================================================================
plt.hist(p213X_radius, bins = 30, range = (0,65), linewidth = 1, edgecolor = 'black',
         color = 'cornflowerblue')
plt.ylim(0, 55)
plt.xlabel('Radius')
plt.ylabel('Frequency')
plt.title('Nevus radii Distribution of P21 mice treated with 3X (75mg/mL) tamoxifen')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P21_3x/04162018_nevusRadiDistribution.pdf')


#import all files for same experiment p50 1X
sample39 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/1339_tamoxifen1x_P50tab.txt')
sample41 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/1341_tamoxifen1x_P50tab.txt')
sample44 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/1344_tamoxifen1x_P50tab.txt')
sample50 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/1350_tamoxifen1x_P50tab.txt')
sample98 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/1398_tamoxifen1x_P50tab.txt')
sample01 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/1401_tamoxifen1x_P50tab.txt')
sample03 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/1403_tamoxifen1x_P50tab.txt')
sample06 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/1406_tamoxifen1x_P50tab.txt')

#index the column of interest which in this case was nevus area
sample39 = sample39.iloc[:,1]
sample41 = sample41.iloc[:,1]
sample44 = sample44.iloc[:,1]
sample50 = sample50.iloc[:,1]
sample98 = sample98.iloc[:,1]
sample01 = sample01.iloc[:,1]
sample03 = sample03.iloc[:,1]
sample06 = sample06.iloc[:,1]

#append all files into one dataframe (df)
p501X = sample39.append([sample41, sample44, sample50, sample98, sample01, sample03, sample06])

#remove NaN
p501X = p501X.dropna()
p501X = p501X[p501X > 0]

#make a histogram
plt.hist(df, bins = 25, linewidth = 1, edgecolor = 'black')
plt.xlabel('Nevus area (sq µm)')
plt.ylabel('Frequency')
plt.title('Nevus area Distribution of P50 mice treated with 1X (25mg/mL) tamoxifen')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/nevusAreaDistribution.pdf')

plt.hist(p501X_radius, bins = 30, range = (0,60), linewidth = 1, edgecolor = 'black',
         color = 'cornflowerblue')
plt.ylim(0, 40)
plt.xlabel('Radius (sq µm)')
plt.ylabel('Frequency')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P50_1x/nevusRadiDistribution.pdf')


#import all files for same experiment p100 1X
sample98 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_1x/1398_tamoxifen1x_P100tab.txt')
sample03 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_1x/1403_tamoxifen1x_P100tab.txt')
sample06 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_1x/1406_tamoxifen1x_P100tab.txt')

#index the column of interest which in this case was nevus area
sample98 = sample98.iloc[:,1]
sample03 = sample03.iloc[:,1]
sample06 = sample06.iloc[:,1]

#append all files into one dataframe (df)
p1001X = sample98.append([sample03, sample06])

#remove NaN
p1001X = p1001X.dropna()
p1001X = p1001X[p1001X > 0]

#make a histogram
plt.hist(df, bins = 10, linewidth = 1, edgecolor = 'black')
plt.xlabel('Nevus area (sq µm)')
plt.ylabel('Frequency')
plt.title('Nevus area Distribution of P100 mice treated with 1X (25mg/mL) tamoxifen')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P100_1x/nevusAreaDistribution.pdf')

plt.hist(p1001X_radius, bins = 15, range = (0,60), linewidth = 1, edgecolor = 'black',
         color = 'cornflowerblue')
plt.ylim(0, 40)
plt.xlabel('Radius (sq µm)')
plt.ylabel('Frequency')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P100_1x/nevusRadiDistribution.pdf')



#import all files for same experiment p100 3X
sample30 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/1130_tamoxifen3x_P100tab.txt')
sample136 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/1136_tamoxifen3x_P100tab.txt')
sample137 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/1137_tamoxifen3x_P100tab.txt')
sample32 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/1232_tamoxifen3x_P100tab.txt')
sample33 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/1233_tamoxifen3x_P100tab.txt')
sample35 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/1235_tamoxifen3x_P100tab.txt')
sample236 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/1236_tamoxifen3x_P100tab.txt')
sample237 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/1237_tamoxifen3x_P100tab.txt')
T296 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P100_tabs/T296_PBS_P100_3x.txt')
T306 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P100_tabs/T306_pbs_P100_3x.txt')
T311 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/MF_P100_tabs/T311_pbs_P100_3x.txt')

#index the column of interest which in this case was nevus area
sample30 = sample30.iloc[:,1]
sample136 = sample136.iloc[:,1]
sample137 = sample137.iloc[:,1]
sample32 = sample32.iloc[:,1]
sample33 = sample33.iloc[:,1]
sample35 = sample35.iloc[:,1]
sample236 = sample236.iloc[:,1]
sample237 = sample237.iloc[:,1]
T296 = T296.iloc[:,1]
T306 = T306.iloc[:,1]
T311 = T311.iloc[:,1]

#append all files into one dataframe (df)
p1003X = sample30.append([sample136, sample137, sample32, sample33, sample35, sample236, sample237])

#remove NaN
p1003X = p1003X.dropna()
p1003X = p1003X[p1003X > 0]

#make a histogram
plt.hist(df, bins = 25, linewidth = 1, edgecolor = 'black')
plt.xlabel('Nevus area (sq µm)')
plt.ylabel('Frequency')
plt.title('Nevus area Distribution of P100 mice treated with 3X (75mg/mL) tamoxifen')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/nevusAreaDistribution.pdf')

plt.hist(p1003X_radius, bins = 15, range = (0,60), linewidth = 1, edgecolor = 'black',
         color = 'cornflowerblue')
plt.ylim(0, 40)
plt.xlabel('Radius (sq µm)')
plt.ylabel('Frequency')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/P100_3x/nevusRadiDistribution.pdf')


p211X_radius = p211X.apply(lambda x: math.sqrt(x/math.pi))
p501X_radius = p501X.apply(lambda x: math.sqrt(x/math.pi))
p1001X_radius = p1001X.apply(lambda x: math.sqrt(x/math.pi))
p1003X_radius = p1003X.apply(lambda x: math.sqrt(x/math.pi))

#import all files for same ATR experiment 3X
sample176 = pd.read_table('/Users/Ruiz/Documents/Ganesan_Lab/Tab_files/ATR/T176_ATRmut3x_MFdepltab.txt')

#index the column of interest which in this case was nevus area
p50 = sample176.iloc[:,8]

#remove NaN
df = p50.dropna()

#make a histogram
plt.hist(df, bins = 10, linewidth = 1, edgecolor = 'black')
plt.xlabel('Nevus area (sq µm)')
plt.ylabel('Frequency')
plt.title('Nevus area Distribution of P50 ATR mice treated with 3X (75mg/mL) tamoxifen')
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/ATR/p50MtnevusAreaDistribution.pdf')

#import all files for same P50 WT control for ATR experiment 3X 
#No nevi


#make graph of average sizes. drop na and drop 0. 
plt.bar(range(3), [p21_mean, p50_mean, p100_mean], width = 0.5, yerr = [p21_err, p50_err, p100_err], color = 'lightsalmon')
plt.ylabel('Average nevi area treated with 3X 4-OHT')
plt.xlabel('Age')
plt.xticks(range(3), ('p21', 'p50', 'p100'))
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/averageNevusSize.pdf')

plt.bar(range(2), [p21_mean, p50_mean], width = 0.6, yerr = [p21_err, p50_err], 
        color = 'cornflowerblue', capsize = 5)
plt.ylabel('Average nevi area treated with 3X 4-OHT (sq µm)')
plt.xlabel('Age')
plt.xticks(range(2), ('p21', 'p50'))
plt.savefig('/Users/ruiz/Documents/Ganesan_Lab/Tab_files/averageNevusSize.pdf')









