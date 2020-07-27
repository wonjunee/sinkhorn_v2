import matplotlib.pyplot as plt
import numpy as np
import dask.dataframe as dd
import pandas as pd
import sys
import cv2



def open_csv(filename,n1,n2):
    A=dd.read_csv(filename,header=None)
    return np.array(A).reshape((n2,n1))

with open("./data/parameters.csv") as F:
    n1=int(F.readline())
    n2=int(F.readline())
    nt=int(F.readline())
    tau=float(F.readline())
    gamma=float(F.readline())
    kappa=float(F.readline())

video_type = 0

if len(sys.argv)>1:
    nt = int(sys.argv[1])

if len(sys.argv)>2:
    video_type = int(sys.argv[2])

# #--------------------------------------------------
# #   Getting Rho Data
# #--------------------------------------------------

def open_and_reshape(filename, n1, n2):
    X = np.fromfile(filename, dtype=np.float64)
    X.shape = (n1,n2)
    return X


#--------------------------------------------------
#   Create animation
#--------------------------------------------------
from matplotlib import animation



x = np.linspace(0, 1, n1)
y = np.linspace(0, 1, n2)
x, y = np.meshgrid(x, y)

if video_type == 0:
    # First set up the figure, the axis, and the plot element we want to animate
    rho=open_and_reshape("./data/mu-{}.csv".format(0),n1,n2)
    fig, ax = plt.subplots(1,2,figsize=(10,5))
    fig.subplots_adjust(bottom=0.05, top=0.95, right=1, left=0)
    cax = ax[1].imshow(rho, cmap='inferno')

    rho_bin = np.zeros_like(rho)
    rho_bin[rho>0] = 1
    cax1= ax[0].imshow(rho_bin, cmap='inferno')
    fig.colorbar(cax)
    ax[0].set_axis_off()
    ax[1].set_axis_off()

    # max_rho = np.max(rho)
    max_rho = np.max(rho)
    # animation function.  This is called sequentially
    def animate(n):

        print("\rProcessing frame %d/%d..." % (n, nt), end='')

        # fig.clear()
        rho=open_and_reshape("./data/mu-{}.csv".format(n),n1,n2)
        # cax.set_array(np.flipud(rho))
        cax.set_array(rho)
        cax.set_clim(np.min(0), np.max(rho))
        rho_bin = np.zeros_like(rho)
        rho_bin[rho>1e-3] = 1
        cax1.set_array(rho_bin)
        cax1.set_clim(0, 1)
        # plt.show()
        
        # cax1.set_clim(np.min(0), np.max(rho))
        plt.suptitle("q={}, tau={}, gamma={}\nn={}, t={:.2f}, sum={:.4f}".format(kappa,tau,gamma,n,n*tau,np.sum(rho)/(n1*n2)),fontsize=6)
        return cax, cax1,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt+1, interval=10, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save("video-with-title.mp4", fps=30, dpi = 200)

elif video_type == 1:
    # First set up the figure, the axis, and the plot element we want to animate
    rho=open_and_reshape("./data/mu-{}.csv".format(0),n1,n2)
    SIZE_RATIO = 100
    fig, ax = plt.subplots(1,1,figsize=(6,6))
    fig.subplots_adjust(bottom=0, top=1, right=1, left=0)
    cax = ax.imshow(rho, cmap='inferno')
    plt.axis('off')

    def animate(n):

        print("\rProcessing frame %d/%d..." % (n, nt), end='')

        # fig.clear()
        rho=open_and_reshape("./data/mu-{}.csv".format(n),n1,n2)
        # cax.set_array(np.flipud(rho))
        # obstacle[((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5) <= 0.2**2)] = np.max(rho)*0.5

        # rho[((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5) <= 0.201**2)] = np.max(rho)*0.5

        cax.set_array(rho)
        cax.set_clim(np.min(0), np.max(rho))

        return cax, 

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt+1, interval=10, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save("video.mp4", fps=30)