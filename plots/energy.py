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
plt.title(fr'Energia con $N = {N}$, $dt = ${dt:.1e}')
textstr = f'$V_0 = ${V0:.1e}\n$a = {a}$'
props = dict(boxstyle='round', facecolor='white')
plt.gcf().text(0.76, 0.78, textstr, fontsize=12, verticalalignment='bottom', bbox=props)
plt.xlabel('t')
plt.ylabel('E')
plt.legend()
plt.tight_layout()

# plt.savefig('report/figures/gauss_energy.png', dpi=500)
plt.show()