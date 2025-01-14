import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from functions import getParam

N, M, L, dx, dt, V0, a = getParam("data/trapped/param.csv")
data = pd.read_csv("data/trapped/energy.csv")

t = data['t']
K = data['K']
V = data['V']
E = data['E']

plt.figure()
plt.plot(t, K, label=r'$\langle K \rangle$')
plt.plot(t, V, label=r'$\langle V \rangle$')
plt.plot(t, E, label=r'$\langle H \rangle$')
plt.title(fr'Energy plot')
plt.xlabel('t')
plt.ylabel('E')
plt.legend()
plt.tight_layout()

# plt.savefig('report/figures/pos.png', dpi=500)
plt.show()