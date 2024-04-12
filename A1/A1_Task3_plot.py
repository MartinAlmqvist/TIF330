import numpy as np
import matplotlib.pyplot as plt

# Load data from the file
data = np.loadtxt('results.txt')

# Extract Ey and Bz from the data
Ey = data[:, 0]
Bz = data[:, 1]

# Create an array for x values
x = np.linspace(-0.1, 0.1, len(Ey))

# Plot Ey and Bz
plt.plot(x, Ey, label='Ey')
plt.plot(x, Bz, label='Bz')
plt.xlabel('x')
plt.ylabel('Field Strength')
plt.title('Electric and Magnetic Field Evolution')
plt.legend()
plt.grid(True)
plt.show()
