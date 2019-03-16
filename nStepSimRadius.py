import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

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

#print(stats.kstest(FinalRadius, 'lognorm'))

plt.hist(FinalRadius, bins=50, linewidth = 1, edgecolor = 'black',
         density=True, color = 'cornflowerblue')
plt.xlabel(r'Radius ($\mu$m)')
plt.ylabel('Frequency')
plt.title('Simulation of %d Steps Transformation with Transformation Probability %.2f' %(n, p))
xmin, xmax = plt.xlim()
lnspc = np.linspace(xmin, xmax, len(FinalRadius))

m1, s1, t1 = stats.lognorm.fit(FinalRadius)
pdf_log = stats.lognorm.pdf(lnspc, m1, s1, t1)
plt.plot(lnspc, pdf_log, label = "Lognorm fit", color = "g", linewidth=3)
plt.legend()
plt.show()
