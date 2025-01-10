import matplotlib.pyplot as plt
import pandas as pd

dataCN = pd.read_csv("data/C-N.csv")
dataEU = pd.read_csv("data/EU.csv")
t = dataCN['t']
norm_sqCN = dataCN['norm_sq']
norm_sqEU = dataEU['norm_sq']

plt.figure()
plt.plot(t, norm_sqCN, label='CN')
plt.plot(t, norm_sqEU, label='EU')
plt.legend()

plt.show()