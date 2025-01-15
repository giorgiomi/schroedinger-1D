import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from functions import getParam

# Load parameters and data for TS simulation
N, M, L, dx, dt, V0, a = getParam("data/test/paramT-S.csv")
data_ts = pd.read_csv("data/test/energy.csv")

# Load parameters and data for CN simulation
N_cn, M_cn, L_cn, dx_cn, dt_cn, V0_cn, a_cn = getParam("data/trapped/param.csv")
data_cn = pd.read_csv("data/trapped/energy.csv")

# Check that the parameters are the same
assert (N, M, L, dx, dt, V0, a) == (N_cn, M_cn, L_cn, dx_cn, dt_cn, V0_cn, a_cn), "Parameters do not match between TS and CN simulations"

# Extract data for TS simulation
t_ts = data_ts['t']
K_ts = data_ts['K']
V_ts = data_ts['V']
E_ts = data_ts['E']

# Extract data for CN simulation
t_cn = data_cn['t']
K_cn = data_cn['K']
V_cn = data_cn['V']
E_cn = data_cn['E']

plt.figure()

# Plot CN simulation data
plt.plot(t_cn, K_cn, label=r'$\langle T \rangle_{CN}$')
plt.plot(t_cn, V_cn, label=r'$\langle V \rangle_{CN}$')
plt.plot(t_cn, E_cn, label=r'$\langle H \rangle_{CN}$')

# Plot TS simulation data
plt.plot(t_ts, K_ts, label=r'$\langle T \rangle_{TS}$', linestyle='--')
plt.plot(t_ts, V_ts, label=r'$\langle V \rangle_{TS}$', linestyle='--')
plt.plot(t_ts, E_ts, label=r'$\langle H \rangle_{TS}$', linestyle='--')

plt.title(fr'Energia con $N = {N}$, $dt = ${dt:.1e}')
textstr = f'$V_0 = ${V0:.1e}\n$a = {a}$'
props = dict(boxstyle='round', facecolor='white')
plt.gcf().text(0.76, 0.78, textstr, fontsize=12, verticalalignment='bottom', bbox=props)
plt.xlabel('t')
plt.ylabel('E')
plt.legend()
plt.tight_layout()


# # Calculate energy differences
# E_diff = E_ts - E_cn

# # Plot energy difference
# plt.figure()
# plt.plot(t_ts, E_diff, label=r'$\Delta E = \langle H \rangle_{TS} - \langle H \rangle_{CN}$')

# plt.title(fr'Differenza di Energia con $N = {N}$, $dt = ${dt:.1e}')
# plt.xlabel('t')
# plt.ylabel(r'$\Delta E$')
# plt.legend()
# plt.tight_layout()

plt.savefig('report/figures/gauss_energy_compare.png', dpi=500)
plt.show()


## PROBABILITIES

# Load parameters and data for TS simulation
N, M, L, dx, dt, V0, a = getParam("data/test/paramT-S.csv")
data_ts = pd.read_csv("data/test/T-S.csv")

# Load parameters and data for CN simulation
N_cn, M_cn, L_cn, dx_cn, dt_cn, V0_cn, a_cn = getParam("data/trapped/param.csv")
data_cn = pd.read_csv("data/trapped/C-N.csv")

# Extract data for TS simulation
t_ts = data_ts['t']
prob_left_ts = data_ts['p_l']
prob_right_ts = data_ts['p_r']

# Extract data for CN simulation
t_cn = data_cn['t']
prob_left_cn = data_cn['p_l']
prob_right_cn = data_cn['p_r']

plt.figure()

# Plot CN simulation data
plt.plot(t_cn, prob_left_cn, label='$P_L^{CN}$')
plt.plot(t_cn, prob_right_cn, label='$P_R^{CN}$')
plt.plot(t_cn, prob_right_cn + prob_left_cn, label='$P_L^{CN} + P_R^{CN}$')

# Plot TS simulation data
plt.plot(t_ts, prob_left_ts, label='$P_L^{TS}$', linestyle='--')
plt.plot(t_ts, prob_right_ts, label='$P_R^{TS}$', linestyle='--')
plt.plot(t_ts, prob_right_ts + prob_left_ts, label='$P_L^{TS} + P_R^{TS}$', linestyle='--')

plt.axhline(0.5, color='tab:red')
props = dict(boxstyle='round', facecolor='white')
textstr = f'$V_0 = ${V0:.1e}\n$a = {a}$'
plt.gcf().text(0.76, 0.78, textstr, fontsize=12, verticalalignment='bottom', bbox=props)
plt.legend()
plt.title(fr'Oscillazioni di probabilità con $N = {N}$, $dt = ${dt:.1e}')
plt.xlabel('t')
plt.ylabel(r'$P(t)$')
plt.yticks([i * 0.1 for i in range(11)])  # Set y ticks to 0.1 intervals
# plt.yscale('log')
plt.tight_layout()

plt.savefig(f'report/figures/gauss_prob_compare.png', dpi=500)
plt.show()