import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# Read the CSV file
data = np.genfromtxt('wave_function_values.csv', delimiter=',')

# Extract time and wave function values
time = data[:, 0]
wave_function = data[:, 1] + 1j * data[:, 2]  # Combine real and imaginary parts

# Compute the Fourier transform of the wave function
fourier_transform = fft(wave_function)

# Compute the frequencies corresponding to the Fourier transform
sampling_frequency = 1 / (time[1] - time[0])
frequencies = fftfreq(len(time), 1 / sampling_frequency)

# Plot the temporal spectrum
plt.figure(figsize=(10, 6))
plt.plot(frequencies, np.abs(fourier_transform))
plt.title('Temporal Spectrum of Psi(0.1, 0)')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.xlim(0, 0.5 * sampling_frequency)  # Plot only positive frequencies
plt.grid(True)
plt.show()
