import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from functions import getParam

N, M, L, dx, dt = getParam("data/free/param.csv")
dataCN = pd.read_csv("data/free/C-N.csv")
# dataEU = pd.read_csv("data/free/EU.csv")

t = dataCN['t']
x_CN = dataCN['x']
x2_CN = dataCN['x2']
# x_EU = dataEU['x']
# x2_EU = dataEU['x2']

plt.figure()
plt.plot(t, x_CN, label=r'$\langle x \rangle$')
# plt.plot(t, x2_CN, label=r'$\langle x^2 \rangle$')
plt.plot(t, np.sqrt(x2_CN - x_CN**2), label=r'$\sigma_x$')
plt.title(fr'position averages N = {N}, L = {L}, dt = {dt}')
plt.xlabel('t')
# plt.ylabel('x')
plt.legend()
plt.tight_layout()

# plt.figure()
# plt.plot(t, x_EU, label='x')
# plt.plot(t, x2_EU, label='x2')
# plt.plot(t, x2_EU - x_EU**2, label='var')
# plt.legend()

plt.savefig('report/figures/pos.png', dpi=500)
plt.show()