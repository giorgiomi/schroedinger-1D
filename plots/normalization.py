import matplotlib.pyplot as plt
import pandas as pd

dataCN = pd.read_csv("data/C-N.csv")
dataEU = pd.read_csv("data/EU.csv")
t = dataCN['t']
norm_sqCN = dataCN['norm_sq']
norm_sqEU = dataEU['norm_sq']

x_CN = dataCN['x']
x2_CN = dataCN['x2']

x_EU = dataEU['x']
x2_EU = dataEU['x2']

plt.figure()
plt.plot(t, norm_sqCN, label='CN')
plt.plot(t, norm_sqEU, label='EU')
plt.legend()
plt.yscale('log')

plt.figure()
plt.plot(t, x_CN, label='x')
plt.plot(t, x2_CN, label='x2')
plt.plot(t, x2_CN - x_CN**2, label='var')
plt.legend()

plt.figure()
plt.plot(t, x_EU, label='x')
plt.plot(t, x2_EU, label='x2')
plt.plot(t, x2_EU - x_EU**2, label='var')
plt.legend()

plt.show()