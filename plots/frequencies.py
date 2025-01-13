import matplotlib.pyplot as plt
import pandas as pd
from functions import getParam

# Load the frequencies data
data = pd.read_csv('data/trapped/frequencies.csv')
N, M, L, dx, dt, V0, a = getParam("data/trapped/param.csv")

# Plot the frequencies as a function of a or V
if 'a' in data:
    plt.figure()
    x = data['a']
    plt.errorbar(x, data['frequency'], yerr=data['error'], marker='o', linestyle='-', capsize=4)
    plt.xlabel('a')
    plt.ylabel(r'$\nu$')
    plt.title(f'Frequenza di oscillazione con $V_0=${V0}')
    plt.tight_layout()

    # plt.savefig('report/figures/freq_vs_a.png')
elif 'V' in data:
    plt.figure()
    x = data['V']
    plt.errorbar(x, data['frequency'], yerr=data['error'], marker='o', linestyle='-', capsize=4)
    plt.xlabel('$V_0$')
    plt.ylabel(r'$\nu$')
    plt.title(f'Frequenza di oscillazione con $a=${a}')
    plt.tight_layout()

    # plt.savefig('report/figures/freq_vs_V0.png', dpi=500)
elif 'N' in data:
    plt.figure()
    x = data['N']
    plt.errorbar(x, data['frequency'], yerr=data['error'], marker='o', linestyle='-', capsize=4)
    plt.xlabel('$N$')
    plt.ylabel(r'$\nu$')
    plt.title(f'Frequenza di oscillazione con $a=${a}')
    plt.tight_layout()

    # plt.savefig('report/figures/freq_vs_N.png', dpi=500)

plt.show()