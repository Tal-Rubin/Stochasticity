from cmath import sqrt
#from re import A
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from scipy.special import comb


import matplotlib.animation as animation
from matplotlib.ticker import FormatStrFormatter


#mpl.style.use("dark_background")
mpl.rcParams['axes.formatter.limits'] = (-2,2)
mpl.rcParams['agg.path.chunksize'] = 10000000  
plt.rcParams["axes.titlesize"] = 14
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 12
if (False):
  plt.rcParams['figure.dpi'] = 300
  plt.rcParams['savefig.dpi'] = 300
  plt.rcParams['axes.titlepad'] = 20

  plt.rcParams["axes.facecolor"] = "black"

protonChtoM = 9.57883322e7
proton_mass = 1.672621911e-27



def poincare_maps(folder_path,num):
  '''  t_num = data[:,1]
  x_t_num = data[:,2]
  y_t_num = data[:,3]
  z_t_num = data[:,4]

  vx_t_num = data[:,5]
  vy_t_num = data[:,6]
  vz_t_num = data[:,7]

  px_t_num = data[:,8]
  py_t_num = data[:,9]
  pz_t_num = data[:,10]'''

  plt.figure()

  for i in range(num):
    header = np.genfromtxt(folder_path+"/ext"+str(i)+".txt", max_rows=1)

    data = np.genfromtxt(folder_path+"/ext"+str(i)+".txt",skip_header=1, skip_footer=1) #id_, t_, x_[0], x_[1], x_[2], v_[0], v_[1], v_[2]

    rho_t_num =   data[:,11]
    theta_t_num = data[:,12]
    Y_t_num =     data[:,13]

  
    plt.scatter(Y_t_num,rho_t_num, s=0.3)
  for i in range(num):
    header = np.genfromtxt(folder_path+"/ext"+str(i)+".txt", max_rows=1)

    data = np.genfromtxt(folder_path+"/ext"+str(i)+".txt",skip_header=1, skip_footer=1) #id_, t_, x_[0], x_[1], x_[2], v_[0], v_[1], v_[2]

    rho_t_num =   data[:,11]
    theta_t_num = data[:,12]
    Y_t_num =     data[:,13]

  
    plt.scatter(Y_t_num[0],rho_t_num[0],color='black', marker='x')

  plt.suptitle("$\\epsilon = {:.2f}, \\alpha$ = {}, $\\omega$ = {}".format(header[1]/header[2],header[1], header[2]))
  plt.xlim(0,2*np.pi)
  plt.ylim(0,2)
  plt.xlabel("Y")
  plt.ylabel("$\\rho$")
  plt.tight_layout()
  plt.savefig("al{}_om{}.pdf".format(header[1], header[2]))

  plt.show()

def plot_xy(folder_path,num):
  header = np.genfromtxt(folder_path+"/ext"+str(num)+".txt", max_rows=1)
  data = np.genfromtxt(folder_path+"/ext"+str(num)+".txt",skip_header=1, skip_footer=1) #id_, t_, x_[0], x_[1], x_[2], v_[0], v_[1], v_[2]

  plt.figure()
  plt.plot(data[:,2],data[:,3])
  ax = plt.gca()
  ax.set_aspect('equal')

  plt.figure()
  plt.subplot(2,1,1)
  plt.plot(data[:,1],data[:,2])
  plt.xlabel("t")
  plt.ylabel("x") 
  plt.subplot(2,1,2)
  plt.plot(data[:,1],data[:,3])
  plt.xlabel("t")
  plt.ylabel("y")


  plt.show()
