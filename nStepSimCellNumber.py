import numpy as np
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

FinalCells = np.array([])
for i in range(5000): # Total number of simulations
    t = 0
    cells = [1]
    for i in range(n):
        cells.append(0)

    while sum(cells)-cells[n] != 0:  # There are still cells not in the final state
        t += 1
        cells = moveCells(cells, p, n)
        if t > CellCycleLimit: # The cell cycle reaches the limit
            break

    FinalCells = np.append(FinalCells, np.sum(cells))
    print(cells)
    print(np.sum(cells))

#print(stats.kstest(FinalRadius, 'lognorm'))

plt.hist(FinalCells, bins=50, linewidth = 1, edgecolor = 'black',
         density=True, color = 'cornflowerblue')
plt.xlabel('Cell Number')
plt.ylabel('Frequency')
plt.title('Simulation of %d Steps Transformation with Transformation Probability %.2f' %(n, p))

plt.show()
