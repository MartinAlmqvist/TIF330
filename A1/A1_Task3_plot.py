import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def set_plot_style():
    sns.set_context("paper", font_scale=2)
    sns.set_style("darkgrid")
    sns.set_palette("deep")
    sns.set(font='sans-serif')

    fig_width = 10
    fig_height = 6

    plt.rcParams['figure.figsize'] = (fig_width, fig_height)
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['lines.linewidth'] = 1
    plt.rcParams['lines.markersize'] = 4
    plt.rcParams['font.size'] = 2.5 * fig_height
    plt.rcParams['xtick.labelsize'] = 2 * fig_height
    plt.rcParams['ytick.labelsize'] = 2 * fig_height
    plt.rcParams['axes.labelsize'] = 2.5 * fig_height
    plt.rcParams['axes.titlesize'] = 3 * fig_height 
    plt.grid(True)

# Set the plot style
set_plot_style()

# Load data 
Ey_data = np.loadtxt('results_Ey.txt')
Bz_data = np.loadtxt('results_Bz.txt')
max_value = np.max(Ey_data)
min_value = np.min(Ey_data)

print(max_value)
print(min_value)
print(Ey_data.shape[1])

# create time and space steps
t = np.arange(0, Ey_data.shape[0]) * 0.01  # Time
x = np.linspace(-1, 1, Ey_data.shape[1])  # Space

# create meshgrid for t and x
T, X = np.meshgrid(t, x)

""" # Plot Ey
plt.pcolormesh(X, T, Ey_data.T, cmap='viridis', shading='auto')
plt.colorbar(label='E-field ($\hat{y}$)')
plt.xlabel('Space ($\hat{x}$)')
plt.ylabel('Time')
plt.title('Change in E-field over Time and Space - Inifinite Space')
plt.tight_layout()
plt.xlim(-0.1, 0.2)
plt.ylim(0,0.5)
plt.savefig('Figures/A1_E_inf_zoom.png')
#plt.savefig('Figures/A1_E_inf.png')
plt.show()

# Plot Bz
plt.pcolormesh(X, T, Bz_data.T, cmap='viridis', shading='auto')
plt.colorbar(label='B-field ($\hat{z}$)')
plt.xlabel('Space ($\hat{x}$)')
plt.ylabel('Time')
plt.title('Change in B-field over Time and Space - Inifinite Space')
plt.tight_layout()
plt.savefig('Figures/A1_B_inf.png')
plt.show()
 """
plt.pcolormesh(X, T, Ey_data.T, cmap='viridis', shading='auto')
plt.colorbar(label='E-field ($\hat{y}$)')
plt.xlabel('Space ($\hat{x}$)')
plt.ylabel('Time')
plt.title('Change in E-field over Time and Space - Ideal Conductor')
plt.tight_layout()
""" plt.xlim(-0.1, 0.2)
plt.ylim(0,0.5) """
#plt.savefig('Figures/A1_E_IC_zoom.png')
plt.savefig('Figures/A1_E_IC_disp.png')
plt.show()

# Plot Bz
plt.pcolormesh(X, T, Bz_data.T, cmap='viridis', shading='auto')
plt.colorbar(label='B-field ($\hat{z}$)')
plt.xlabel('Space ($\hat{x}$)')
plt.ylabel('Time')
plt.title('Change in B-field over Time and Space - Ideal Conductor')
plt.tight_layout()
plt.savefig('Figures/A1_B_IC_disp.png')
plt.show()  