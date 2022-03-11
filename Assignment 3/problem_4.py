import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

np.random.seed(2021)
X = np.random.rand(6)
X = np.sort(X)
Y = np.random.rand(6)

cubic_spline = CubicSpline(X, Y, bc_type='natural', extrapolate=True)
plt.scatter(X, Y,marker='o', label='data', color='red')
plt.plot(np.linspace(0, 1-1e-20, 200), cubic_spline(np.linspace(0, 1-1e-20, 200)),  label="interpolation", color='blue')
plt.legend(loc='upper right')
plt.title('Natural Cubic Interpolation',color='#641E16')
plt.xlabel('X-values',color='#7E5109')
plt.ylabel('Y-values',color='#7E5109')
plt.show()
#plt.savefig('problem_4c.png')

