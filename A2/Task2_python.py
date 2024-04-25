import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# Load data from CSV file
data = np.loadtxt('wave_function_values.csv', delimiter=',')

# Extract time and wave function values
time = data[:, 0]
wave_function = data[:, 1]

# Compute FFT of the wave function
fft_values = fft(wave_function)
freqs = fftfreq(len(time), d=time[1] - time[0])

# Plot temporal spectrum
plt.figure()
plt.plot(freqs, np.abs(fft_values))
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.title('Temporal Spectrum of Wave Function at (0.1, 0)')
plt.grid(True)
plt.show()
