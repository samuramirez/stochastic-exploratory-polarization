
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import zipfile
import shutil

import scipy.io
import os
import os.path
from scipy.optimize import curve_fit
from functions_patch import get_centroid
from functions_patch import sum_squared_displacements_sphere
import pickle
plt.ion()


matplotlib.rcParams.update({'font.size': 7})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1


folder='D:/dynamic_polarity_data/fig_3D_transition_quant/'
subfolder = 'newreac_gef100_k4a'
parameters=['0p001','0p002','0p005','0p01','0p02','0p05','0p1','1']
filenames = ['traj_3d_newreac_gef100_k4a_'+parameters[i]+'_' for i in range(len(parameters))]
skipframes=[2+3,2+3,2+3,2+3,2+3,2+3,2+3,2+3]
R=4.5135/2 

#extract files
if not os.path.isdir(folder+subfolder):
    print ('unzipping files\n')
    with zipfile.ZipFile(folder+subfolder + '.zip',"r") as zip_ref:
        zip_ref.extractall(folder)
    
    print ('done unzipping files\n')

        

labels= [str(parameters[i]) for i in range(len(parameters))]
#colors=['blue','red','green','black']
#colors = [ cm.jet(0.2), cm.jet(0.5) , cm.jet(0.7) ,cm.jet(0.9) ]

colors = [cm.hot(i) for i in np.linspace(0.5, 0, len(parameters))]
 
Nsims=30
T=1*3600
sampling= 1 #1min intervals for centroid calculation
sampling_sim = 60 #s sampling time in simulations
interval=T/sampling_sim #nframes


#MSD
#max step in a trajectory to be considered part of the same path to estimate 
maxstep=6
sum_norm_disps=np.zeros((len(parameters),interval)) 
sum_squared_disps=np.zeros((len(parameters),interval)) 
sum_squared_disps_squared=np.zeros((len(parameters),interval))
n_squared_disps=np.zeros((len(parameters),interval))


for p,filename in enumerate(filenames):
    
    #nframes=T/sampling_sim - skipframes[p]
    #samples=int(nframes/sampling)
    for sim in range(Nsims):
        print(sim)
        name =folder + subfolder +'/' + filename + '{:02d}'.format(sim+1) + '.txt'
      
        if os.path.isfile(name):
            #print(name)
            TS_centroid = np.genfromtxt(name,skip_header=skipframes[p])                
 
            sum_norm_disp, sum_squared_disp, sum_squared_disp_squared, n_squared_disp = sum_squared_displacements_sphere(TS_centroid, maxstep, R)       
            sum_norm_disps[p,0:len(sum_norm_disp)] +=  sum_norm_disp   
            sum_squared_disps[p,0:len(sum_squared_disp)] +=  sum_squared_disp  
            sum_squared_disps_squared[p,0:len(sum_squared_disp_squared)] +=  sum_squared_disp_squared  
            n_squared_disps[p,0:len(n_squared_disp)] +=  n_squared_disp
            
          
        else:
            print "file " + name + " not found"
           
data = {'ssd':sum_squared_disps, 'nsq':n_squared_disps, 'params':parameters,
         'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
        'mxstp':maxstep}

with open(folder+subfolder+'.pkl', 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open('filename.pickle', 'rb') as handle:
 #   b = pickle.load(handle)

#erase folder 
#shutil.rmtree(folder+subfolder)


