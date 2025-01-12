import matplotlib.pyplot as plt
import pandas as pd
from functions import getParam
from scipy.signal import find_peaks

N, M, L, dx, dt, V0, a = getParam("data/trapped/param.csv")
dataCN = pd.read_csv("data/trapped/C-N.csv")

t = dataCN['t']
prob_left = dataCN['p_l']
prob_right = dataCN['p_r']
# print(prob_left.to_string())

# Find intersection points between prob_left and the horizontal line 0.5
intersection_indices = []
for i in range(1, len(prob_left)):
    if (prob_left[i-1] - 0.5) * (prob_left[i] - 0.5) < 0:
        intersection_indices.append(i)
# print("Indices of intersection points with 0.5 in prob_left:", intersection_indices)

# Transform the indices into times
intersection_times = [t[i] for i in intersection_indices]

# Calculate the time differences between consecutive intersection times
time_differences = [intersection_times[i] - intersection_times[i-1] for i in range(1, len(intersection_times))]
# Keep only the time differences greater than the first one/2
time_differences = [diff for diff in time_differences if diff > time_differences[0] / 2]

# Calculate the average frequency
if time_differences:
    average_period = sum(time_differences) / len(time_differences)
    average_frequency = 1 / average_period
else:
    average_frequency = 0

print(f"Average Frequency: {average_frequency:.3f}")

plt.figure()

plt.plot(t, prob_left, label='Left')
plt.plot(t, prob_right, label='Right')
plt.plot(t, prob_right + prob_left)

if time_differences:
    for i, (time, diff) in enumerate(zip(intersection_times, time_differences)):
        if diff > time_differences[0] / 2:
            plt.annotate('', xy=(time, 0.5), xytext=(time, 0.8),
                         arrowprops=dict(facecolor='black', edgecolor='none', shrink=0.05, width=1, headwidth=5),
                         label='Intersections' if i == 0 else "")
    plt.annotate('', xy=(intersection_times[i] + time_differences[-1], 0.5), xytext=(intersection_times[i] + time_differences[-1], 0.8),
                 arrowprops=dict(facecolor='black', edgecolor='none', shrink=0.05, width=1, headwidth=5))

plt.axhline(0.5, color='tab:red')
# Add a box with the parameters V0 and a
props = dict(boxstyle='round', facecolor='white')
textstr = f'$V_0 = {V0}$\n$a = {a}$\n'
textstr += fr'$\langle \nu \rangle = {average_frequency:.3f}$'
plt.gcf().text(0.3, 0.18, textstr, fontsize=12, verticalalignment='bottom', bbox=props)
plt.legend()
plt.title(fr'Oscillation N = {N}, L = {L}, dt = {dt}')
plt.xlabel('t')
plt.ylabel(r'$P(t)$')
# plt.yscale('log')
plt.tight_layout()


# plt.savefig('report/figures/norm.png', dpi=500)
plt.show()