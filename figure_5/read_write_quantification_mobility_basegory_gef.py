
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


folder='./'


subfolder = 'basegory_newreac0_gef'

parameters=[450,500,600,700]
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[15,15,15,15] #


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
Nsims=50
T=1*3600
sampling= 1 #1min intervals for centroid calculation
sampling_sim = 60 #s sampling time in simulations
interval=T/sampling_sim #nframes

centerwrt=1 #center with respect to frame
Nspecies=1
L=8.0


#MSD
#max step in a trajectory to be considered part of the same path to estimate 
maxstep=6
sum_norm_disps=np.zeros((len(parameters),interval)) 
sum_squared_disps=np.zeros((len(parameters),interval)) 
sum_squared_disps_squared=np.zeros((len(parameters),interval))
n_squared_disps=np.zeros((len(parameters),interval))
stepshistograms=[]
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
    samples=int(nframes/sampling)
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
            
            TS_tot42T = TS_cdc42T_gef42
            
            #CENTER
            #get coordinates of max value in the first frame
            mat_first = np.reshape(TS_tot42T[centerwrt,:],(N,N))
            im,jm = np.unravel_index(mat_first.argmax(), mat_first.shape)
                        
            TS_centered=np.zeros((samples,N*N))
            for i in range(samples):
                #set as matrix
                mat=np.reshape(TS_tot42T[i*sampling,:],(N,N))
                #center
                centered=mat[np.ix_((np.arange(N) + im - int(N/2)) % N , (np.arange(N) + jm - int(N/2)) % N)]
                #store as vector
                TS_centered[i,:]=centered.flat
              
            #Compute centroid trajectory 
            TS_centroid=[]
            for i in range(samples):
                mat = np.reshape(TS_centered[i,:],(N,N))
                #if there are proteins
                if np.sum(mat>0) > 0:
                    [xcent, ycent] = get_centroid(mat)
                    TS_centroid.append([xcent, ycent])
                else:
                    TS_centroid.append([np.nan, np.nan])
        
            TS_centroid=np.asarray(TS_centroid)   
            
            #get steps
            dx = (TS_centroid[1:,0] - TS_centroid[0:-1,0])*L/N
            dy = (TS_centroid[1:,1] - TS_centroid[0:-1,1])*L/N
            #boundary conditions
            dx[dx < -L/2] = L + dx[dx < -L/2] 
            dx[dx > L/2] = -L + dx[dx > L/2]
            dy[dy < -L/2] = L + dy[dy < -L/2] 
            dy[dy > L/2] = -L + dy[dy > L/2]
            
            steps = np.sqrt(dx**2 + dy**2 )
            #if np.sum(steps>1.5) > 0:
            #    print (sim,steps)
            stepsizes=np.concatenate((stepsizes,steps))
            #vector steps
            vec_step = TS_centroid[1:]-TS_centroid[0:-1]
            #Angle [-pi,pi]
            theta=np.arctan2(vec_step[:,1],vec_step[:,0])
            #set nan when the step size is too small
            #theta[(stepsize <= np.sqrt(2)*pixel_size )] = np.nan
            #theta[(stepsize <= pixel_size )] = np.nan
            #thetas=np.concatenate((thetas,theta))  
            
            #turning angle    
            turnp=theta[1:]-theta[0:-1]
            #turning angle [-pi,pi]
            turn=np.arctan2(np.sin(turnp),np.cos(turnp))
            turns=np.concatenate((turns,turn))
                
            sum_norm_disp, sum_squared_disp, sum_squared_disp_squared, n_squared_disp = sum_squared_displacements(TS_centroid*h, maxstep)       
            sum_norm_disps[p,0:len(sum_norm_disp)] +=  sum_norm_disp   
            sum_squared_disps[p,0:len(sum_squared_disp)] +=  sum_squared_disp  
            sum_squared_disps_squared[p,0:len(sum_squared_disp_squared)] +=  sum_squared_disp_squared  
            n_squared_disps[p,0:len(n_squared_disp)] +=  n_squared_disp
            
            #Tot 42T
            #sum over grid element            
            tot42tvst = np.sum(TS_tot42T[:,:],1)
            #sum over time            
            tot42t+=np.sum(tot42tvst)            
            tot42t2+=np.sum(tot42tvst**2)
            
          
            
            ntot42t+=Ntimes            
        else:
            print "file " + name + " not found"
           
    stepshistograms.append(stepsizes)
    meantot42t.append(1.0*tot42t/ntot42t)
    meantot42t2.append(1.0*tot42t2/ntot42t)
  


data = {'ssd':sum_squared_disps, 'nsq':n_squared_disps, 'params':parameters,
        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
         'mxstp':maxstep}

with open(folder+subfolder+'.pkl', 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open('filename.pickle', 'rb') as handle:
 #   b = pickle.load(handle)

#erase folder 
#shutil.rmtree(folder+subfolder)


