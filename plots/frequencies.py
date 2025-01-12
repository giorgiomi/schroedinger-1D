import matplotlib.pyplot as plt
import pandas as pd

# Load the frequencies data
data = pd.read_csv('data/trapped/frequencies.csv')

# Plot the frequencies as a function of a
plt.figure(figsize=(10, 6))
plt.errorbar(data['a'], data['frequency'], yerr=data['error'], marker='o', linestyle='-', capsize=5)
plt.xlabel('a')
plt.ylabel('Frequency')
plt.title('Frequency as a function of a')
plt.grid(True)
# plt.savefig('plots/frequency_vs_a.png')
plt.show()