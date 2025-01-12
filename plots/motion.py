import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.animation import FuncAnimation
from functions import getParam

data = pd.read_csv("data/trapped/C-N.csv")
N, M, L, dx, dt = getParam("data/trapped/param.csv")

# Extract time and psi data
time = data.iloc[:, 0].values
psi_real = data.iloc[:, 1:-3:2].values
psi_imag = data.iloc[:, 2:-3:2].values

# Insert values of 0 at the beginning and at the end of each array
psi_real = np.hstack([np.zeros((psi_real.shape[0], 1)), psi_real, np.zeros((psi_real.shape[0], 1))])
psi_imag = np.hstack([np.zeros((psi_imag.shape[0], 1)), psi_imag, np.zeros((psi_imag.shape[0], 1))])
x = np.arange(psi_real.shape[1])*dx - L

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

ani = FuncAnimation(fig, update, frames=len(time), blit=True, interval=100)
plt.xlabel('x')
plt.ylabel('Psi')
plt.title(f'motion plot N = {N}, L = {L}, dt = {dt:.2e}')
plt.ylim(-10, 10) 
plt.show()