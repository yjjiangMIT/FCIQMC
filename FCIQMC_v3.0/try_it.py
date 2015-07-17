import matplotlib
matplotlib.rcParams['backend'] = 'Qt4Agg'
import matplotlib.pyplot as plt
# import pylab as plt
import numpy as np

# Sample data
sidex = [1,2,3,4]
sidey = [1,2,3,4]
X,Y = np.meshgrid(sidex,sidey)
Z = np.array([[1,2,3,0],[4,5,6,0],[7,8,9,0],[0,0,0,0]])

# Plot the density map using nearest-neighbor interpolation
plt.figure(10)
plt.ion()
plt.pcolormesh(X,Y,Z)
plt.colorbar()
plt.xlim([1,4])
plt.ylim([4,1])
plt.show()