import numpy as np
import matplotlib.pyplot as plt

# Load Ey data from the file
Ey_data = np.loadtxt('results_Bz.txt')

# Create a meshgrid for time and space
t = np.arange(0, Ey_data.shape[0]) * 0.01  # Time
x = np.linspace(-0.1, 0.1, Ey_data.shape[1])  # Space

# Create meshgrid for plotting
T, X = np.meshgrid(t, x)

# Plot meshgrid
plt.figure(figsize=(10, 6))
plt.pcolormesh(T, X, Ey_data.T, cmap='viridis')  # Transpose Ey_data
plt.colorbar(label='Ey')
plt.xlabel('Time')
plt.ylabel('Space')
plt.title('Meshgrid for Ey over Time and Space')
plt.show()
