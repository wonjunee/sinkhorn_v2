import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import csv


def open_file(filename):
    arr = []
    with open(filename) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            arr.append(list(map(float,row)))
    return np.array(arr)

xarr = [4, 8, 16, 32, 64, 128, 256, 512, 1024]

def plot_in_plot(filename, ax, label, color):
    arr = open_file(filename)
    arr = arr[1:]
    arr_mean = np.mean(arr, axis = 1)
    arr_max = np.max(arr, axis = 1)
    arr_min = np.min(arr, axis = 1)
    arr_std_up   = arr_mean + np.std(arr, axis = 1)
    arr_std_down = arr_mean - np.std(arr, axis = 1)

    print(filename)
    print(arr_mean)


    ax.loglog(xarr, arr_mean, 'o-', c=color, label=label)
    # plt.fill_between(xarr, arr_min, arr_max,
    #                  facecolor="orange", # The fill color
    #                  color='blue',       # The outline color
    #                  alpha=0.15)          # Transparency of the fill
    ax.fill_between(xarr, arr_min, arr_max,
                     facecolor="orange", # The fill color
                     color=color,       # The outline color
                     alpha=0.2)          # Transparency of the fill

    xarr_np = np.array(xarr)
    # N = 100
    # ax.loglog(xarr, np.power(N, -2.0/xarr_np))



fig, ax = plt.subplots(1,1,figsize=(10,5))

plot_in_plot("error_N_100.dat", ax, "$N=100$", 'r')
plot_in_plot("error_N_500.dat", ax, "$N=500$", "g")
plot_in_plot("error_N_2000.dat", ax, "$N=2000$", "b")


theoretical_rate = np.array([0.1 * np.power(100, -2.0/n) for n in xarr])
ax.loglog(xarr, theoretical_rate, '--', label="theoretical rate", c="magenta")

ax.set_xlabel("Dimension")
ax.set_ylabel("Error")
ax.set_xlim([25,500])
ax.set_xticks(xarr)
ax.set_xticklabels(xarr)
ax.set_yticks([10,1,0.1,0.01,0.001,0.0005])
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("dim_N_100.png")

