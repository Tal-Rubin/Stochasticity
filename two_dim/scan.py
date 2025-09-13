import numpy as np
import os

import post

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

  Lz = 0

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
  for j in range(0):
    initial_conditions.append({'rho':rho0+j*rho_step,'Y':Y0+j*Y_step, 'id_num':j+i+1})


  for ic in initial_conditions:
    rho = ic['rho']
    Y = ic['Y']
    id_num = ic['id_num']
    commandStr = "./poincare_map {} {} {} {} {} {} {} {} {} {} {}".format(len(ions),alpha, dimless_wave_freq, B0, theta_sample, species_index, tstep_fac, X, Y, rho, theta0)

    for ion in ions:
      commandStr += " {} {} {}".format(ion['Z'], ion['m'], ion['n'])
    commandStr += " > "+folder_path+"/ext" + str(id_num) + ".txt"

    print(commandStr)

    os.system(commandStr)

#  header = np.genfromtxt("ext.txt", max_rows=1)
#  data = np.genfromtxt("ext.txt",skip_header=2, skip_footer=1) #id_, t_, x_[0], x_[1], x_[2], v_[0], v_[1], v_[2]
  
  post.poincare_maps(folder_path,len(initial_conditions))