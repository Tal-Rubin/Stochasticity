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


def integrate(integrand,r):
    res = np.zeros(integrand.shape)
    for i in range(res.shape[0]):
        res[i] = np.trapz(integrand[0:i],r[0:i])
    return res


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def plot(folder_path,num):
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

  plt.suptitle("$\\alpha$ = {}, $\\omega$ = {}".format(header[1], header[2]))
  plt.xlim(0,2*np.pi)
  plt.ylim(0,2)
  plt.xlabel("Y")
  plt.ylabel("$\\rho$")
  plt.tight_layout()
  plt.savefig("al{}_om{}.pdf".format(header[1], header[2]))

  plt.show()


if __name__=='__main__':
  tstep_fac = 4000
  theta_sample = -np.pi/2

  alpha = 0.09
  dimless_wave_freq = -0.19

  folder_path = 'alpah{}_om{}'.format(alpha, dimless_wave_freq)
  if not os.path.exists(folder_path):
      os.makedirs(folder_path)

  B0 = 10

  ions = []    
  ions.append({'Z':1,'m':1,'n':.5})
  ions.append({'Z':5,'m':11,'n':.1})

  species_index = 0


  initial_conditions = []

  Y0 = 0.94
  X = 0
  rho0 = 0.065
  theta0 = theta_sample
  rho_step = 0.06
  Y_step = 0.045
  for i in range(30):
    initial_conditions.append({'rho':rho0+i*rho_step,'Y':Y0+i*Y_step, 'id_num':i})
  
  i=29
  rho0 = 0.065
  Y0 = 1.8
  rho_step = 0.06
  #Y_step = 0.03
  for j in range(15):
    initial_conditions.append({'rho':rho0+j*rho_step,'Y':Y0+j*Y_step, 'id_num':j+i+1})


  for ic in initial_conditions:
    rho = ic['rho']
    Y = ic['Y']
    id_num = ic['id_num']
    commandStr = "./plasmaSlabSim {} {} {} {} {} {} {} {} {} {} {}".format(len(ions),alpha, dimless_wave_freq, B0, theta_sample, species_index, tstep_fac, X, Y, rho, theta0)

    for ion in ions:
      commandStr += " {} {} {}".format(ion['Z'], ion['m'], ion['n'])
    commandStr += " > "+folder_path+"/ext" + str(id_num) + ".txt"

    print(commandStr)

    os.system(commandStr)

#  header = np.genfromtxt("ext.txt", max_rows=1)
#  data = np.genfromtxt("ext.txt",skip_header=2, skip_footer=1) #id_, t_, x_[0], x_[1], x_[2], v_[0], v_[1], v_[2]
  
  plot(folder_path,len(initial_conditions))