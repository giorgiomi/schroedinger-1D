import matplotlib.pyplot as plt
import pandas as pd
from functions import getParam

N, M, L, dx, dt = getParam("data/trapped/param.csv")
dataCN = pd.read_csv("data/trapped/C-N.csv")

t = dataCN['t']
prob_left = dataCN['p_l']
prob_right = dataCN['p_r']

plt.figure()
plt.plot(t, prob_left, label='Left')
plt.plot(t, prob_right, label='Right')
plt.legend()
plt.title(fr'Oscillation N = {N}, L = {L}, dt = {dt}')
plt.xlabel('t')
plt.ylabel(r'$P(t)$')
# plt.yscale('log')
plt.tight_layout()


# plt.savefig('report/figures/norm.png', dpi=500)
plt.show()