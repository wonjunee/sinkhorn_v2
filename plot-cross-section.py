import matplotlib.pyplot as plt
import matplotlib
# from mayavi import mlab
import numpy as np
import dask.dataframe as dd
import pandas as pd
import sys
from skimage import data, color
from skimage.transform import rescale, resize, downscale_local_mean

from matplotlib import animation
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import os
from matplotlib import cm

matplotlib.use('Agg')

def open_csv(filename,n1,n2):
    A=dd.read_csv(filename,header=None)
    return np.array(A).reshape((n2,n1))

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


# class Exact:
#     def __init__(self,n1,n2,tau,m):
#         self.n1=n1
#         self.n2=n2
#         self.tau=tau
#         self.m=m
#         self.mprime = self.m/(self.m-1)
#         self.d= 2

#         self.alpha = self.d/(self.d*(self.m-1) + 2)
#         self.k = (self.m-1)*self.alpha/(2.0*self.m*self.d)
#         self.C = ( self.mprime * self.k / np.pi ) ** (1/self.mprime)


#         # set xi0 as a free parameter
#         self.xi0 =  0.2
#         self.s = ( self.xi0 * self.xi0 * self.k / self.C ) ** (1/self.alpha)

#         # set s as a free parameter
#         self.s   = 1e-4
#         self.xi0 = ( self.C * self.s ** self.alpha / self.k ) ** 0.5

#         # set initial height as a free parameter
#         self.h0 = 5
#         self.s  = (self.mprime / (np.pi*self.h0**m)) ** (1.0/(-self.alpha*self.m - self.alpha - 2))
#         self.C  = (self.h0 * self.s**(-self.alpha)) ** (self.m - 1)

#         print("self.s = ", self.s)
        

#     def F(self,x):
#         return np.maximum(self.C - self.k * x * x, 0) ** (1.0/(self.m - 1))

#     def U(self,x,y,t):
#         r = np.sqrt(x*x+y*y)
#         return ((t + self.s))**(-self.alpha) * self.F( r * ((t + self.s))**(-self.alpha/self.d) )


class Exact:
    def __init__(self,n1,n2,tau,m,M):
        self.n1=n1
        self.n2=n2
        self.tau=tau
        self.m=m
        self.mprime = self.m/(self.m-1)
        self.d= 2
        self.M=M

        self.alpha = self.d/(self.d*(self.m-1) + 2)
        self.k = (self.m-1)*self.alpha/(2.0*self.m*self.d)
        self.C = (self.M * self.mprime * self.k / np.pi ) ** (1/self.mprime)

        # set initial height as a free parameter
        self.h0 = 15
        self.s  = (self.h0 * self.C**(-1.0/(self.m-1))) ** (-1.0/self.alpha)

        print("self.s = ", self.s)
        

    def F(self,x):
        return np.maximum(self.C - self.k * x * x, 0) ** (1.0/(self.m - 1))

    def U(self,x,y,t):
        r = np.sqrt(x*x+y*y)
        return ((t + self.s))**(-self.alpha) * self.F( r * ((t + self.s))**(-self.alpha/self.d) )


class Norm:
    def __init__(self, n1,n2):
        self.n1=n1
        self.n2=n2

    def calc_L1(self,u):
        return np.sum(np.abs(u))/(self.n1*self.n2)

    def calc_L2(self,u):
        return np.sqrt(np.sum(np.power(u,2))/(self.n1*self.n2))

#--------------------------------------------------
#   Create animation
#--------------------------------------------------

# animation function.  This is called sequentially

def make_movie_3d(vmax, frames_dir, num_frames, n1, n2, M):


    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    dpi = n1 / 5.0

    # ax.set_axis_off()
    # fig.subplots_adjust(bottom=0, top=0.8, right=1, left=0)

    ex = Exact(n1,n2,tau,m,M)
    xx = np.linspace(-0.5,0.5,n1,endpoint=False)
    yy = np.linspace(-0.5,0.5,n2,endpoint=False)
    xx = xx+0.5/n1
    yy = yy+0.5/n1
    xx, yy = np.meshgrid(xx, yy)

    
    norm = Norm(n1,n2)

    cumulative_norm = 0

    for k in range(num_frames+1):

        print("\rProcessing frame %d/%d..." % (k, num_frames), end='')

        # mu = open_and_reshape(filename+".data", n1, n2)
        mu = open_and_reshape("./data/mu-{}.csv".format(k),n1,n2)

        if k==0:
            vmax=np.max(mu)

        plt.cla()
        cax = ax.plot(np.linspace(0,1,n1), mu[n2//2,:],label="Computed")

        if k>=0:

            zz = ex.U(xx,yy,k*tau)
            ax.plot(np.linspace(0,1,n1), zz[n2//2,:],'--',label="Exact")

        ax.set_ylim([-0.4,vmax])
        # ax.set_axis_off()
        plt.legend(fontsize=11)

        diff = mu-zz
        L1_norm = norm.calc_L1(diff)
        L2_norm = norm.calc_L2(diff)

        cumulative_norm += L1_norm

        # plt.suptitle("\ntau={}, gamma={}, t={:.3f}".format(tau,gamma,k*tau,L1_norm,L2_norm),fontsize=8)
        plt.suptitle("\ntau={}, gamma={}, t={:.5f}\n$\\||\\rho\\||_{{L^1}}$ = {:.5e}\t$\\||\\rho\\||_{{L^2}}$ = {:.5e}".format(tau,gamma,k*tau,L1_norm,L2_norm),fontsize=11)

        filename='figures/figure-{:0>3d}.eps'.format(k)
        plt.savefig(filename, dpi=dpi)

    print("")
    print("cumulative_norm   : ", cumulative_norm*tau)
    print("cumulative_norm/T : ", cumulative_norm*tau/(tau*nt) )
    print("done")   
    print("")


# ================================================================================================

if __name__ == "__main__":

    with open("./data/parameters.csv") as F:
        n1=int(F.readline())
        n2=int(F.readline())
        nt=int(F.readline())
        tau=float(F.readline())
        gamma=float(F.readline())
        m=float(F.readline())
        M=float(F.readline())


    if len(sys.argv)==2:
        nt = int(sys.argv[1])

    print(n1,n2,nt,tau,gamma,m,M)

    max_rho = 5
    os.system("rm ./figures/figures-*")


    make_movie_3d(max_rho, "", nt, n1, n2, M)

    os.system("ffmpeg -r 10 -f image2 -s 1280x1280 -i ./figures/figure-%03d.png -loglevel quiet -vcodec libx264 -crf 25  -pix_fmt yuv420p cross.mp4")