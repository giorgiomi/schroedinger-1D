import matplotlib.pyplot as plt
import pandas as pd
from functions import getParam
import csv
import warnings
import sys
warnings.filterwarnings("ignore")


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
        if not intersection_indices or (i - intersection_indices[-1]) > intersection_indices[0] / 2:
            intersection_indices.append(i)
# print("Indices of intersection points with 0.5 in prob_left:", intersection_indices)

# Transform the indices into times
intersection_times = [t[i] for i in intersection_indices]

# Calculate the time differences between consecutive intersection times
time_differences = [intersection_times[i] - intersection_times[i-1] for i in range(1, len(intersection_times))]
# Keep only the time differences greater than the first one/2
time_differences = [diff for diff in time_differences if diff > time_differences[0] / 2]
k = len(time_differences)

# Calculate the average frequency and its error
if time_differences:
    T_avg = sum(time_differences) / k
    f_avg = 1 / T_avg
    # Calculate the standard deviation of the periods
    T_std = (sum((diff - T_avg) ** 2 for diff in time_differences) / k) ** 0.5
    # Calculate the error on the average frequency
    df = T_std / ((k ** 0.5) * (T_avg ** 2))
else:
    f_avg = 0
    df = 0

# Save the results to a CSV file
with open('data/trapped/frequencies.csv', 'a', newline='') as csvfile:
    writer = csv.writer(csvfile)
    if sys.argv[2] == 'a':
        writer.writerow([a, f_avg, df])
    elif sys.argv[2] == 'V0':
        writer.writerow([V0, f_avg, df])
    elif sys.argv[2] == 'N':
        writer.writerow([N, f_avg, df])

if sys.argv[1] != 'y':
    exit()

plt.figure()

plt.plot(t, prob_left, label='$P_L$')
plt.plot(t, prob_right, label='$P_R$')
plt.plot(t, prob_right + prob_left, label='$P_L + P_R$')

# if time_differences:
#     for i, (time, diff) in enumerate(zip(intersection_times, time_differences)):
#         if diff > time_differences[0] / 2:
#             plt.annotate('', xy=(time, 0.5), xytext=(time, 0.8),
#                          arrowprops=dict(facecolor='black', edgecolor='none', shrink=0.05, width=1, headwidth=5),
#                          label='Intersections' if i == 0 else "")
#     plt.annotate('', xy=(intersection_times[i] + time_differences[-1], 0.5), xytext=(intersection_times[i] + time_differences[-1], 0.8),
#                  arrowprops=dict(facecolor='black', edgecolor='none', shrink=0.05, width=1, headwidth=5))

plt.axhline(0.5, color='tab:red')
# Add a box with the parameters V0 and a
props = dict(boxstyle='round', facecolor='white')
textstr = f'$V_0 = ${V0:.1e}\n$a = {a}$\n'
textstr += fr'$\langle \nu \rangle = {f_avg:.1f}\pm {df:.1f}$'
plt.gcf().text(0.76, 0.74, textstr, fontsize=12, verticalalignment='bottom', bbox=props)
plt.legend()
plt.title(fr'Oscillazioni di probabilità con $N = {N}$, $dt = ${dt:.1e}')
plt.xlabel('t')
plt.ylabel(r'$P(t)$')
plt.yticks([i * 0.1 for i in range(11)])  # Set y ticks to 0.1 intervals
# plt.yscale('log')
plt.tight_layout()


# plt.savefig(f'report/figures/prob_null.png', dpi=500)
plt.show()