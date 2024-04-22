import numpy as np
import matplotlib.pyplot as plt

def xfun(x):
    return 2*x / (4*(1/4 + x**2))

x = np.linspace(-6, 6, 1000)
y = xfun(x)

max_y = np.max(y)
max_x = x[np.argmax(y)]

plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.plot(x, y, label=r'$\frac{b}{2(\frac{1}{4} + b^2)}$')
plt.scatter(max_x, max_y, color='red', label=r'$\sigma_{\mathrm{max}}$: (1/2, 1/2)')
plt.xlabel('b',fontsize=16)
plt.ylabel('$\sigma(b)$',fontsize=16)
plt.title('Sigma as a function of b with $a^2$= 1/4',fontsize=16)
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
plt.legend(fontsize=14)
plt.show()
