
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
from functions_patch import sum_squared_displacements
import pickle
plt.ion()


matplotlib.rcParams.update({'font.size': 7})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1


folder='D:/dynamic_polarity_data/basegory/'


subfolder = 'basegory_k1a_10x_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k2a_10_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k2b_0p0625x_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k3_10_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k4a_10_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k5a_100x_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k7_50_gef'
parameters=[15,25,50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k1b_0p0625x_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k5b_0p0625x_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing

subfolder = 'basegory_k4b_0p0625x_gef'
parameters=[50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10] # skipping a lot of frames for high gef due to multiple patches competing






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
 
N1=80
Nsims=3
T=600
sampling_sim = 30 #s sampling time in simulations
interval=T/sampling_sim #nframes
centerwrt=1 #center with respect to frame
Nspecies=2
L=8.0


meantot42t=[]
meantot42t2=[] 
meangef42t=[]
meangef42t2=[] 

for p,filename in enumerate(filenames):
    stepsizes=np.empty(0)
    turns=np.empty(0)
    TS_cdc42T_gef42=0 #empty cdc42t distribution to avoid reusing when a .dat file is not present
    tot42t=0
    tot42t2=0
    gef42t=0
    gef42t2=0
    ntot42t=0
    nframes=T/sampling_sim - skipframes[p]
    #samples=int(nframes/sampling)
    for sim in range(Nsims):
        print(sim)
        name =folder + subfolder +'/' + filename + str(sim+1) + '.dat'
      
        if os.path.isfile(name):
            #print(name)
            TS_cdc42T_gef42 = np.genfromtxt(name,skip_header=skipframes[p]*Nspecies)
            Nrows=np.shape(TS_cdc42T_gef42)[0]    
            Ntimes=Nrows/Nspecies
            if Ntimes != nframes:
                print "Incorrect number of frames"
                print(str(Ntimes), str(nframes), str(sim))
                
            N=int(np.sqrt(np.shape(TS_cdc42T_gef42)[1]))
            if N1 != N:
                print "Incorrect input N"
            h=L/N
            
            #Get totalCdc42T = Cdc42T+GEF42 
            TS_tot42T = TS_cdc42T_gef42[0::2,:] + TS_cdc42T_gef42[1::2,:]
            TS_gef42T = TS_cdc42T_gef42[1::2,:]
            
            
            #Tot 42T
            #sum over grid element            
            tot42tvst = np.sum(TS_tot42T[:,:],1)
            #sum over time            
            tot42t+=np.sum(tot42tvst)            
            tot42t2+=np.sum(tot42tvst**2)
            
            #GEF42
            gef42tvst = np.sum(TS_gef42T[:,:],1)
            #sum over time            
            gef42t+=np.sum(gef42tvst)            
            gef42t2+=np.sum(gef42tvst**2)
            
            ntot42t+=Ntimes            
        else:
            print "file " + name + " not found"
           
    meantot42t.append(1.0*tot42t/ntot42t)
    meantot42t2.append(1.0*tot42t2/ntot42t)
    meangef42t.append(1.0*gef42t/ntot42t)
    meangef42t2.append(1.0*gef42t2/ntot42t)



data = { 'params':parameters,
        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
        'mgef42t': meangef42t, 'mgef42t2': meangef42t2,
        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
        }

with open(folder+subfolder+'.pkl', 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open('filename.pickle', 'rb') as handle:
 #   b = pickle.load(handle)

#erase folder 
#shutil.rmtree(folder+subfolder)


