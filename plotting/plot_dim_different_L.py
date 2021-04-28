import numpy as np
import matplotlib
import pandas as pd
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import csv


def open_file(filename):
    return np.array(pd.read_csv(filename, header=None))

def plot_in_plot(filename, ax, label, color):
    arr = open_file(filename)
    arr_mean = np.mean(arr, axis = 1)
    arr_max = np.max(arr, axis = 1)
    arr_min = np.min(arr, axis = 1)
    arr_std_up   = arr_mean + np.std(arr, axis = 1)
    arr_std_down = arr_mean - np.std(arr, axis = 1)

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

fig, ax = plt.subplots(1,1,figsize=(10,5))
s = "$O(d^{\\lambda}/\\sqrt{N})$"
xarr = [2**(n+1) for n in range(7)]

N=500
val = plot_in_plot("L_error_N_{}.dat".format(N), ax, "original", 'tab:blue')

val = plot_in_plot("L1_error_N_{}.dat".format(N), ax, "$L(\\nu) - L(\\rho)$ L: Gaussian", 'tab:orange')

val = plot_in_plot("L2_error_N_{}.dat".format(N), ax, "$L(\\nu) - L(\\rho)$ L: $1/x$", 'tab:green')

val = plot_in_plot("L3_error_N_{}.dat".format(N), ax, "$L(\\nu) - L(\\rho)$ L: $1/x^2$", 'tab:red')

ax.set_xlabel("Dimension")
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
plt.savefig("dim_N_500_different_L.png")

