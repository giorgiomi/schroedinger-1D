import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.animation import FuncAnimation

data = pd.read_csv("data/C-N.csv")
# Extract time and psi data
time = data.iloc[:, 0].values
psi_real = data.iloc[:, 1:-1:2].values
psi_imag = data.iloc[:, 2:-1:2].values
x = np.arange(psi_real.shape[1])

fig, ax = plt.subplots()
line_real, = ax.plot(x, psi_real[0, :], label='Real part')
line_imag, = ax.plot(x, psi_imag[0, :], label='Imaginary part')
line_norm, = ax.plot(x, psi_real[0, :]**2 + psi_imag[0, :]**2, label='Normalization')
ax.legend()

def update(frame):
    line_real.set_ydata(psi_real[frame, :])
    line_imag.set_ydata(psi_imag[frame, :])
    line_norm.set_ydata(psi_real[frame, :]**2 + psi_imag[frame, :]**2)
    return line_real, line_imag, line_norm

ani = FuncAnimation(fig, update, frames=len(time), blit=True, interval=100)
plt.xlabel('Position')
plt.ylabel('Psi')
plt.title('Motion plot of Psi')
plt.ylim(-1, 1) 
plt.show()