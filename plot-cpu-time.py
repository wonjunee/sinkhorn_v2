import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def open_file(filename):
    arr = []
    with open(filename) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            arr.append(list(map(float,row)))
    return np.array(arr)

xarr = [25, 50, 100, 250, 500]
# yarr = [0.000887, 0.001321, 0.002352, 0.00516, 0.010408]
yarr = [0.000853, 0.001384, 0.002826, 0.004933, 0.009699]

xarr = [25,100,500]
# xarr_power = [1, 2, 2.3, 2.7, 3, 3.3, 3.7, 4]
xarr_power = [1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.3, 3.7]
xarr = [10**i for i in xarr_power]
import csv

arr = open_file('time_check_0.005000.dat')
arr_mean = np.mean(arr, axis = 1)
arr_max = np.max(arr, axis = 1)
arr_min = np.min(arr, axis = 1)
arr_std_up   = arr_mean + np.std(arr, axis = 1)
arr_std_down = arr_mean - np.std(arr, axis = 1)

arr_H = open_file('time_check_0.010000.dat')
arr_H_mean = np.mean(arr_H, axis = 1)
arr_H_max = np.max(arr_H, axis = 1)
arr_H_min = np.min(arr_H, axis = 1)
arr_H_std_up   = arr_H_mean + np.std(arr_H, axis = 1)
arr_H_std_down = arr_H_mean - np.std(arr_H, axis = 1)

arr_H2 = open_file('time_check_0.100000.dat')
arr_H2_mean = np.mean(arr_H2, axis = 1)
arr_H2_max = np.max(arr_H2, axis = 1)
arr_H2_min = np.min(arr_H2, axis = 1)
arr_H2_std_up   = arr_H2_mean + np.std(arr_H2, axis = 1)
arr_H2_std_down = arr_H2_mean - np.std(arr_H2, axis = 1)

print(arr)
fig, ax = plt.subplots(1,1)

ax.loglog(xarr, arr_mean, 'o-', c="blue", label="$\\lambda=0.005$")
# plt.fill_between(xarr, arr_min, arr_max,
#                  facecolor="orange", # The fill color
#                  color='blue',       # The outline color
#                  alpha=0.15)          # Transparency of the fill
plt.fill_between(xarr, arr_std_down, arr_std_up,
                 facecolor="orange", # The fill color
                 color='blue',       # The outline color
                 alpha=0.2)          # Transparency of the fill


ax.loglog(xarr, arr_H_mean, 'o-', c="red", label="$\\lambda=0.01$")
# plt.fill_between(xarr, arr_H_min, arr_H_max,
#                  facecolor="orange", # The fill color
#                  color='red',       # The outline color
#                  alpha=0.15)          # Transparency of the fill
plt.fill_between(xarr, arr_H_std_down, arr_H_std_up,
                 facecolor="orange", # The fill color
                 color='red',       # The outline color
                 alpha=0.2)          # Transparency of the fill

ax.loglog(xarr, arr_H2_mean, 'o-', c="orange", label="$\\lambda=0.1$")
# plt.fill_between(xarr, arr_H2_min, arr_H2_max,
#                  facecolor="orange", # The fill color
#                  color='red',       # The outline color
#                  alpha=0.15)          # Transparency of the fill
plt.fill_between(xarr, arr_H2_std_down, arr_H2_std_up,
                 facecolor="orange", # The fill color
                 color='red',       # The outline color
                 alpha=0.2)          # Transparency of the fill


theoretical_rate = np.array([0.15 * np.power(n, -0.4) for n in xarr])
ax.loglog(xarr, theoretical_rate, '--', label="theoretical rate", c="magenta")

ax.set_xlabel("Number of points")
ax.set_ylabel("Error")
ax.set_xlim([25,500])
ax.set_xticks(xarr)
ax.set_xticklabels(["10^{}".format(i) for i in xarr_power])
# ax.set_yticks([0.1,0.01,0.001,0.0005])
ax.grid()
ax.legend()
plt.savefig("excution-time-H-hat-and-normal.png")
plt.show()


# fig, ax = plt.subplots(1,1)
# ax.loglog(xarr, yarr, 'o-')
# ax.set_xlabel("Dimension")
# ax.set_ylabel("Execution time in seconds")
# ax.set_xlim([25,500])
# ax.set_xticks(xarr)
# ax.set_xticklabels(xarr)
# ax.set_yticks([0.1,0.01,0.001,0.0005])
# ax.grid()
# # plt.savefig("excution-time.png")
# plt.show()