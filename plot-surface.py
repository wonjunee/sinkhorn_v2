import argparse

import matplotlib.pyplot as plt
# from mayavi import mlab
import numpy as np
import dask.dataframe as dd
import pandas as pd
import sys
from skimage import data, color
from skimage.transform import rescale, resize, downscale_local_mean

def open_csv(filename,n1,n2):
    A=dd.read_csv(filename,header=None)
    return np.array(A).reshape((n2,n1))

with open("./data/parameters.csv") as F:
    n1=int(F.readline())
    n2=int(F.readline())
    nt=int(F.readline())
    tau=float(F.readline())
    gamma=float(F.readline())


parser = argparse.ArgumentParser(description='Short sample app')

parser.add_argument('-n', action='store', dest='nt',
                    help='Store value of nt')

parser.add_argument('-v', action='store', dest='video_type',
                    help='Store value of video type. 0: with title 1: without title, default: with title')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

if(results.nt != None):
    nt = int(results.nt)
video_type = int(results.video_type)

print('nt           =', nt)
print('video_type   =', video_type)

print(n1,n2,nt,tau,gamma)

def plot_mu():
    
    mu = open_csv("./data/mu.csv",n1,n2)
    nu = open_csv("./data/nu.csv",n1,n2)
    phi = open_csv("./data/phi.csv",n1,n2)
    psi = open_csv("./data/psi.csv",n1,n2)
    push_mu = open_csv("./data/push_mu.csv",n1,n2)

    fig, ax = plt.subplots(1,3)
    ax[0].imshow(mu,cmap='inferno')
    ax[0].set_axis_off()
    ax[0].set_title("mu")
    ax[1].imshow(nu,cmap='inferno')
    ax[1].set_axis_off()
    ax[1].set_title("nu")
    ax[2].imshow(push_mu,cmap='inferno')
    ax[2].set_axis_off()
    ax[2].set_title("push_mu")
    plt.show()

    fig, ax = plt.subplots(1,2)
    ax[0].imshow(phi,cmap='inferno')
    ax[0].set_axis_off()
    ax[0].set_title("phi")
    ax[1].imshow(psi,cmap='inferno')
    ax[1].set_axis_off()
    ax[1].set_title("psi")
    plt.show()


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
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

X = np.linspace(0, 1, n1)
Y = np.linspace(0, 1, n2)
X, Y = np.meshgrid(X, Y)

# surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# Customize the z axis.

max_rho = 5

import os
from matplotlib import cm
os.system("rm ./figures/surface-*")

# animation function.  This is called sequentially



def make_movie_3d(vmax, frames_dir, num_frames, n1, n2):

    SIZE = 256

    fig = plt.figure(figsize=(6,4))  # 16/9 ratio
    ax = fig.add_subplot(111, projection='3d')
    dpi = 960 / 8.0

    ax.set_axis_off()
    fig.subplots_adjust(bottom=-0.15, top=1, right=1.35, left=-0.40)

    x = np.linspace(0, 1, SIZE)
    y = np.linspace(0, 1, SIZE)
    x, y = np.meshgrid(x, y)

    for k in range(num_frames+1):

        # ax.view_init(15, (-k*0.5) % 360)
        # ax.view_init(20, (360-0.5*k)%360)
        ax.view_init(70, 80)

        print("\rProcessing frame %d/%d..." % (k, num_frames), end='')

        # mu = open_and_reshape(filename+".data", n1, n2)
        mu = open_and_reshape("./data/mu-{}.csv".format(k),n1,n2)

        if k==0:
            vmax=np.max(mu)*0.8


        # mu = resize(mu, (SIZE, SIZE),
        #                anti_aliasing=True)
        mu = resize(mu, (SIZE, SIZE))
        plt.cla()
        # cax = ax.plot_surface(x, y, mu, cmap='inferno', vmin=0.0, vmax=vmax)
        cax = ax.plot_surface(x, y, mu, rstride=6, cstride=6, cmap='inferno',
                        vmin=0, vmax=vmax)
        ax.set_zlim(0.0, vmax)

        ax.set_axis_off()

        if(video_type == 0):
            plt.suptitle("tau={}, gamma={}\nn={}, t={:.2f}".format(tau,gamma,k,k*tau),fontsize=6)

        filename='figures/surface-{:0>3d}.png'.format(k)
        plt.savefig(filename, dpi=dpi)

    print("done")   

make_movie_3d(max_rho, "", nt, n1, n2)

os.system("ffmpeg -loglevel panic -r 10 -f image2 -s 1280x1280 -i ./figures/surface-%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test2.mp4")