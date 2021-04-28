import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import csv
import pandas as pd

dat_folder = "data"

def open_file(filename):
    return np.array(pd.read_csv(filename, header=None))

def plot_in_plot(filename, ax, label, color):
    arr = open_file(filename)
    arr_mean = np.mean(arr, axis = 1)
    arr_max = np.max(arr, axis = 1)
    arr_min = np.min(arr, axis = 1)
    arr_std_up   = arr_mean + np.std(arr, axis = 1)
    arr_std_down = arr_mean - np.std(arr, axis = 1)

    print(filename)
    print(arr_mean)


    ax.plot(xarr, arr_mean, 'o-', c=color, label=label)
    # plt.fill_between(xarr, arr_min, arr_max,
    #                  facecolor="orange", # The fill color
    #                  color='blue',       # The outline color
    #                  alpha=0.15)          # Transparency of the fill
    ax.fill_between(xarr, arr_std_down, arr_std_up,
                     facecolor="orange", # The fill color
                     color=color,       # The outline color
                     alpha=0.2)          # Transparency of the fill
    return arr_mean[0]

def calculate_theoretical_rate(DIM, xarr, val):
    # return np.array([val / np.power(xarr[0], -2.0/DIM) * np.power(n, -2.0/DIM) for n in xarr])

    k = 0
    eps = 0.5
    xarr = np.array(xarr)
    return 0.01*np.exp(k/eps)/np.sqrt(xarr) * (1.0 + np.power(eps, -DIM/10) )

fig, ax = plt.subplots(1,1,figsize=(10,5))

xarr = np.array([100*2**i for i in range(5)])
s = "$O(1/\\sqrt{n})$"

# DIM = 2

# val = plot_in_plot("error_DIM_{}.dat".format(DIM), ax, "d={}".format(DIM), 'r')
# theoretical_rate = calculate_theoretical_rate(DIM, xarr, val)
# ax.plot(xarr, theoretical_rate, '--', label="{} d={}".format(s, DIM))

# DIM = 8
# val = plot_in_plot("error_DIM_{}.dat".format(DIM), ax, "d={}".format(DIM), 'r')
# theoretical_rate = calculate_theoretical_rate(DIM, xarr, val)
# s = "something"
# ax.plot(xarr, theoretical_rate, '--', label="{} d={}".format(s, DIM))

# DIM = 16
# val = plot_in_plot("error_DIM_{}.dat".format(DIM), ax, "d={}".format(DIM), 'g')
# theoretical_rate = calculate_theoretical_rate(DIM, xarr, val)
# s = "something"
# ax.plot(xarr, theoretical_rate, '--', label="{} d={}".format(s, DIM))

DIM = 32



val = plot_in_plot("{}/original_error_DIM_{}.dat".format(dat_folder, DIM), ax, "original".format(DIM), 'tab:green')
val = plot_in_plot("{}/inverse_error_DIM_{}.dat".format(dat_folder, DIM), ax, "inverse".format(DIM), 'tab:blue')
val = plot_in_plot("{}/gaussian_error_DIM_{}.dat".format(dat_folder, DIM), ax, "gaussian 0.1".format(DIM), 'tab:red')
val = plot_in_plot("{}/gaussian_0.5_error_DIM_{}.dat".format(dat_folder, DIM), ax, "gaussian 0.5".format(DIM), 'tab:orange')

theoretical_rate = calculate_theoretical_rate(DIM, xarr, val)
ax.plot(xarr, theoretical_rate, '--', label="{} d={}".format(s, DIM))

# DIM = 32
# val = plot_in_plot("{}/L2_error_DIM_{}.dat".format(dat_folder, DIM), ax, "$L(\\nu) - L(\\rho)$ L: Gaussian", 'tab:green')

# DIM = 32
# val = plot_in_plot("{}/L1_error_DIM_{}.dat".format(dat_folder, DIM), ax, "$L(\\nu) - L(\\rho)$ L: $1/x$", 'tab:orange')

# DIM = 32
# val = plot_in_plot("{}/L3_error_DIM_{}.dat".format(dat_folder, DIM), ax, "$L(\\nu) - L(\\rho)$ L: $1/x^2$", 'tab:red')

ax.set_xlabel("N")
ax.set_ylabel("Error")
# ax.set_xlim([25,500])
ax.set_xticks(xarr)
ax.set_xticklabels(xarr)
# ax.set_yticks([10,1,0.1,0.01,0.001,0.0005])
ax.grid()
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
plt.tight_layout()
plt.savefig("dimplot.png")

