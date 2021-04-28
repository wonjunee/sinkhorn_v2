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

def plot_in_plot(filename, ax, label, color):
    arr = open_file(filename)
    arr_mean = np.mean(arr, axis = 1)
    arr_max = np.max(arr, axis = 1)
    arr_min = np.min(arr, axis = 1)
    arr_std_up   = arr_mean + np.std(arr, axis = 1)
    arr_std_down = arr_mean - np.std(arr, axis = 1)

    print(filename)
    print(arr_mean)

    l = len(xarr)

    ax.plot(xarr, arr_mean[:l], 'o-', c=color, label=label)
    # plt.fill_between(xarr, arr_min, arr_max,
    #                  facecolor="orange", # The fill color
    #                  color='blue',       # The outline color
    #                  alpha=0.15)          # Transparency of the fill
    ax.fill_between(xarr, arr_std_down[:l], arr_std_up[:l],
                     facecolor="orange", # The fill color
                     color=color,       # The outline color
                     alpha=0.2)          # Transparency of the fill
    return arr_mean[0]
eps = 0.9

def calculate_theoretical_rate(N, val, xarr):
    # return np.array([val / np.power(N, -2.0/xarr[0]) * np.power(N, -2.0/n) for n in xarr])
    k = 1
    xarr = np.array(xarr)
    print(xarr)
    # return 0.008*np.exp(k/eps)/np.sqrt(N) * (0.0 + np.power(eps, -xarr/50) )
    return 0.04/np.sqrt(N)*np.power(xarr,0.5)

fig, ax = plt.subplots(1,1,figsize=(10,5))

xarr = [2**(n+1) for n in range(7)]

N = 100
val = plot_in_plot("error_N_{}.dat".format(N), ax, "N={}".format(N), 'r')
theoretical_rate = calculate_theoretical_rate(N, val, xarr)
s = "$O(d^{\\lambda}/\\sqrt{N})$"
ax.plot(xarr, theoretical_rate, '--', label="{} N={}".format(s,N))

N = 300
val = plot_in_plot("error_N_{}.dat".format(N), ax, "N={}".format(N), 'g')
theoretical_rate = calculate_theoretical_rate(N, val, xarr)
s = "$O(d^{\\lambda}/\\sqrt{N})$"
ax.plot(xarr, theoretical_rate, '--', label="{} N={}".format(s,N))

# N = 1000
# val = plot_in_plot("error_N_{}.dat".format(N), ax, "N={}".format(N), 'b')
# theoretical_rate = calculate_theoretical_rate(N, val, xarr)
# s = "$O(d^{\\lambda}/\\sqrt{N})$"
# ax.plot(xarr, theoretical_rate, '--', label="{} N={}".format(s,N))

N = 900
val = plot_in_plot("error_N_{}.dat".format(N), ax, "N={}".format(N), 'y')
theoretical_rate = calculate_theoretical_rate(N, val, xarr)
s = "$O(d^{\\lambda}/\\sqrt{N})$"
ax.plot(xarr, theoretical_rate, '--', label="{} N={}".format(s,N))

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
plt.savefig("sinkhorn_x_axis_dim_lambda_{}.png".format(eps))

