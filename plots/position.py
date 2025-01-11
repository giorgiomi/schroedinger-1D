import matplotlib.pyplot as plt
import pandas as pd
from functions import getParam

N, M, L, dx, dt = getParam("data/param.csv")
dataCN = pd.read_csv("data/C-N.csv")
# dataEU = pd.read_csv("data/EU.csv")

t = dataCN['t']
x_CN = dataCN['x']
x2_CN = dataCN['x2']
# x_EU = dataEU['x']
# x2_EU = dataEU['x2']

plt.figure()
plt.plot(t, x_CN, label=r'$\langle x \rangle$')
plt.plot(t, x2_CN, label=r'$\langle x^2 \rangle$')
plt.plot(t, x2_CN - x_CN**2, label=r'$\langle x^2 \rangle$ - $\langle x \rangle^2$')
plt.title(fr'position averages N = {N}, L = {L}, dt = {dt}')
plt.xlabel('t')
plt.ylabel('x')
plt.legend()

# plt.figure()
# plt.plot(t, x_EU, label='x')
# plt.plot(t, x2_EU, label='x2')
# plt.plot(t, x2_EU - x_EU**2, label='var')
# plt.legend()

plt.show()