import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("data/psi.csv")
t = data['t']
norm_sq = data['norm_sq']
re = data['re']
im = data['im']

plt.figure()
plt.plot(t, norm_sq)

plt.figure()
plt.plot(t, re, label=r'$Re(\Psi)$')
plt.plot(t, im, label=r'$Im(\Psi)$')
plt.plot(t, re**2 + im**2, label=r'$|\Psi|^2$')
plt.legend()

plt.show()