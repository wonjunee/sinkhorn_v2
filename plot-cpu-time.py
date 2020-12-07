import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


xarr = [25, 50, 100, 250, 500]
yarr = [0.000887, 0.001321, 0.002352, 0.00516, 0.010408]

fig, ax = plt.subplots(1,1)
ax.loglog(xarr, yarr, 'o-')
ax.set_xlabel("Dimension")
ax.set_ylabel("Execution time in seconds")
ax.set_xlim([25,500])
ax.set_xticks(xarr)
ax.set_xticklabels(xarr)
ax.set_yticks([1,0.1,0.01,0.001])
ax.grid()
# plt.savefig("excution-time.png")
plt.show()