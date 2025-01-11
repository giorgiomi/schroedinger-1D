import matplotlib.pyplot as plt
import pandas as pd
from functions import getParam

N, M, L, dx, dt = getParam("data/param.csv")
dataCN = pd.read_csv("data/C-N.csv")
dataEU = pd.read_csv("data/EU.csv")

t = dataCN['t']
norm_sqCN = dataCN['norm_sq']
norm_sqEU = dataEU['norm_sq']

plt.figure()
plt.plot(t, norm_sqEU, label='Euler (explicit)')
plt.plot(t, norm_sqCN, label='Crank-Nicolson')
plt.legend()
plt.title(fr'$\psi$ normalization N = {N}, L = {L}, dt = {dt}')
plt.xlabel('t')
plt.ylabel(r'$\int|\psi(x)|^2dx$')
# plt.yscale('log')

plt.show()