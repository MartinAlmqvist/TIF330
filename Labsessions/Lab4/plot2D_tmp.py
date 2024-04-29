import matplotlib.pyplot as plt
import numpy as np
import struct
f = open('plot2D_tmp.dat', 'rb')
Ny = struct.unpack('I', f.read(4))[0]
Nx = struct.unpack('I', f.read(4))[0]
v = np.empty(shape=(Ny, Nx))
for y in range(Ny):
    for x in range(Nx):
        v[y][x] = struct.unpack('d', f.read(8))[0]
fig, ax = plt.subplots()
vAbsMax = np.maximum(abs(np.amin(v)), abs(np.amax(v)))
plot = ax.imshow(v, interpolation='none', extent=[0,1,0,1], origin='lower', cmap='RdBu', vmin = -vAbsMax, vmax = vAbsMax)
ax.set_aspect('equal')
fig.colorbar(plot, ax=ax, location='right')
plt.gca().set_axis_off()
plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
plt.margins(0,0)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.savefig('T9.png', bbox_inches = 'tight', pad_inches = 0)
