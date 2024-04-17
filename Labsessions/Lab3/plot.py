import numpy as np
import matplotlib.pyplot as plt

y = np.loadtxt('output.txt')
N = 100
x = np.linspace(-1,1, N)

plt.plot(x[1:], y[1:])
plt.show()