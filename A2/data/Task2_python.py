import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# Load the wave function data from the CSV file
data = np.genfromtxt('wave_function_values_t5000_3512.csv', delimiter=',')



# Extract values
time = data[:, 0]
real_part = data[:, 1]
imaginary_part = data[:, 2]

# sin = np.sin(time)
# dt = time[1]-time[0]
# N = len(time)
# yf = fft(sin)
# xf = fftfreq(N,dt)[:N//2]*(2*np.pi)
# print(dt)

sin = (real_part + 1j * imaginary_part)
dt = time[1]-time[0]
N = len(time)
yf = fft(sin)
xf = fftfreq(N,dt)[:N//2]*(2*np.pi)

plt.figure(figsize=(10, 5))
plt.plot(xf, 2.0/N *abs(yf[0:N//2]))
plt.axvline(x = 3.232, color = "r", linestyle = "--", label = "Peak at: 3.232")
plt.xlabel('Frequency ')
plt.ylabel('Absolute magnitude |$\Psi$|')
plt.title('Frequency Spectrum of Wave Function at $\psi(0.1,0)$')
plt.grid()
plt.legend()
plt.show()



# positive_indices = frq >= 0
# positive_frequency = frq[positive_indices]
# temporal_spectrum = np.fft.fft(abs(sin))
# positive_temporal_spectrum = temporal_spectrum[positive_indices]

# plt.plot(time,sin)
# plt.plot(positive_frequency,np.abs(positive_temporal_spectrum))

# magnitude = np.abs(real_part + 1j * imaginary_part)
# dt # plt.figure(figsize=(10, 5))
# #plt.plot(frequency, np.abs(magnitude))
# plt.plot(positive_frequency, np3.232.abs(positive_temporal_spectrum))
# plt.xlabel('Frequency ')
# plt.ylabel('Absolute magnitude [|$\Psi$|]')
# plt.title('Positive Temporal Spectrum of Wave Function at $\psi(0.1,0)$')poral_spectrum = np.fft.fft(magnitude)
# frequency = np.fft.fftfreq(len(time), dt)

# # Filter out positive frequencies and corresponding magnitudes
# positive_indices = frequency >= 0
# positive_frequency = frequency[positive_indices]
# positive_temporal_spectrum = temporal_spectrum[positive_indices]

# # Plot the spectrum
# plt.figure(figsize=(10, 5))
# #plt.plot(frequency, np.abs(magnitude))
# plt.plot(positive_frequency, np.abs(positive_temporal_spectrum))
# plt.xlabel('Frequency ')
# plt.ylabel('Absolute magnitude [|$\Psi$|]')
# plt.title('Positive Temporal Spectrum of Wave Function at $\psi(0.1,0)$')
# plt.grid(True)
#plt.xlim(0, 2)


# # Plot the Wave
# plt.figure(figsize=(10, 5))
# plt.plot(time, magnitude)
# #plt.plot(positive_frequency, np.abs(positive_temporal_spectrum))
# plt.xlabel('Time [s]')
# plt.ylabel('Magnitude [|$\psi(0.1,0)|$]')
# plt.title('Wavefunction over time at $\psi(0.1,0)$')
# plt.grid(True)
# #plt.savefig("Task3_XXYY_wave")
# plt.show()
