import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.animation import FuncAnimation
from functions import getParam
from mpl_toolkits.mplot3d import Axes3D

# data = pd.read_csv("data/trapped/C-N.csv")
# N, M, L, dx, dt, V0, a = getParam("data/trapped/param.csv")
data = pd.read_csv("data/trapped/T-S.csv")
N, M, L, dx, dt, V0, a = getParam("data/trapped/paramT-S.csv")

# Extract time and psi data
time = data.iloc[:, 0].values
psi_real = data.iloc[:, 1:-3:2].values
psi_imag = data.iloc[:, 2:-3:2].values

# Insert values of 0 at the beginning and at the end of each array
psi_real = np.hstack([np.zeros((psi_real.shape[0], 1)), psi_real, np.zeros((psi_real.shape[0], 1))])
psi_imag = np.hstack([np.zeros((psi_imag.shape[0], 1)), psi_imag, np.zeros((psi_imag.shape[0], 1))])
x = np.arange(psi_real.shape[1])*dx - L


## 2D PLOT
fig, ax = plt.subplots()
line_real, = ax.plot(x, psi_real[0, :], label=r'$Re(\psi)$')
line_imag, = ax.plot(x, psi_imag[0, :], label=r'$Im(\psi)$')
line_norm, = ax.plot(x, psi_real[0, :]**2 + psi_imag[0, :]**2, label=r'$|\psi|^2$')
ax.legend()

# Add a text annotation for the time
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def update(frame):
    line_real.set_ydata(psi_real[frame, :])
    line_imag.set_ydata(psi_imag[frame, :])
    line_norm.set_ydata(psi_real[frame, :]**2 + psi_imag[frame, :]**2)
    time_text.set_text(f'Time = {time[frame]:.4f}')
    return line_real, line_imag, line_norm, time_text

ani = FuncAnimation(fig, update, frames=len(time), blit=True, interval=17)

plt.xlabel('x')
plt.ylabel('Psi')
plt.title(f'N = {N}, dt = {dt:.2e}, V_0 = {V0}, a = {a}')
plt.ylim(-2, 2) 
plt.show()

# ani.save(f'data/animations/motion_{N}_{dt:.3f}_{V0}_{a}.mp4', writer='pillow', dpi=300)

exit()

## 3D PLOT
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def update(frame):
    ax.clear()
    ax.plot(x, psi_real[frame, :], zs=psi_imag[frame, :], zdir='z', label=r'$Re(\psi)$')
    ax.plot(x, psi_imag[frame, :], zs=psi_real[frame, :], zdir='y', label=r'$Im(\psi)$')
    ax.set_xlabel('x')
    ax.set_ylabel(r'Re($\psi$)')
    ax.set_zlabel(r'Im($\psi$)')
    ax.set_title(f'motion plot N = {N}, L = {L}, dt = {dt:.2e}')
    ax.set_zlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.legend()
    return ax

ani = FuncAnimation(fig, update, frames=len(time), blit=False, interval=50)
plt.show()
