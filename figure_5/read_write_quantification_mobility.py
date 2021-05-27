
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


#folder = './newreac_mm_gef/'
folder='./'


#subfolder = 'basegory_newreac_k4a_0_gef'
#parameters=[100,150, 200, 300,500]

#subfolder = 'basegory_newreac_gef'
#parameters=[10,15,25,50,75,100,200,300,500]

#subfolder = 'basegory_newreac_k4a_0_k3_0_gef'
#parameters=[100,200, 300, 400,500]

#subfolder = 'basegory_newreac_k4a_0_k3_0_k2a_0p01_gef'
#parameters=[100,200, 300, 400]

#subfolder = 'basegory_k7_k4a_0_gef'
#parameters=[100,150, 200, 300,500,1000]

#subfolder = 'basegory_k7_5_gef'
#parameters=[100,300,500,1000]

#subfolder = 'basegory_k7_5_k3_1_gef'
#parameters=[25,50]
#filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]


#subfolder = 'basegory_newreac_k3_0_gef'
#parameters=[15,25,50,100,200,300,500]
#filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
#skipframes=25 #need to skip competition times

subfolder = 'basegory_k7_1_gef'
parameters=[200,300,500] #100 gef does not polarize with k7 =1 as seen in folder basegory_gef100_k7
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=10 #


subfolder = 'basegory/basegory_k7_10_gef_long'
#parameters=[10,20,50,100,200,500] #20 gef does not polarize with k7 =10 as seen in folder basegory_gef100_k7
parameters=[25,50,100,200,500]
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,25] #

subfolder = 'basegory/basegory_k7_20_gef_long'
#parameters=[10,20,50,100,200,500] #20 gef does not polarize with k7 =10 as seen in folder basegory_gef100_k7
parameters=[25,50,100,200,500]
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,10,20] #

subfolder = 'basegory/basegory_k7_50_gef_long'
parameters=[15,25,50,100,200,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,15,30] #


subfolder = 'basegory_k7_100_gef'
parameters=[10,20,50,100,200,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,20,30] #

subfolder = 'basegory/basegory_k2b_0p0625x_k5a_10x_gef_long'
parameters=[15,25,50,100,300,500] #
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10,10,10] #

subfolder = 'basegory_newreac_gef50_k2b'
parameters=[0.0394, 0.0788, 0.1575, 0.315, 1,3,10]
filenames = ['dynN80k2b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[15,15,15,10,5,5,5]

subfolder = 'basegory_newreac_gef100_k2b'
parameters=[0.0788, 0.1575, 0.315, 1,3,10,30]
filenames = ['dynN80k2b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[20,20,20,10,5,5,5]

subfolder = 'basegory_newreac_gef300_k2b'
parameters=[0.0788, 0.1575, 0.315, 1,3,10,30]
filenames = ['dynN80k2b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[20,20,20,10,5,5,5]

subfolder = 'basegory_newreac_k4a_0_gef'
parameters=[100,150, 200, 300,500]
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5]

subfolder = 'basegory_newreac_k3_1_gef'
parameters=[15,50,100,300,500]
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,15,15,25]

subfolder = 'basegory_newreac_k4a_0_k3_0_gef'
parameters=[100,200,300,400,500]
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5]

subfolder = 'basegory_newreac_gef300_k4a_0_k3'
parameters=[0.005,0.01,0.02,0.04,0.07,0.1,1]
filenames = ['dynN80k3'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5,5,5]

subfolder = 'newreac_controls/gef300_k4b'
parameters=[1,3,30,100]
filenames = ['dynN80k4b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5]

subfolder = 'newreac_controls/gef400_k4a_0_k4b'
parameters=[1,3,30,100]
filenames = ['dynN80k4b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5]


subfolder = 'basegory_newreac_gef200_k4a0_k2b'
parameters=[0.05, 0.1, 0.3,1,3,10]
filenames = ['dynN80k2b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[20,20,10,5,5,5]

subfolder = 'basegory_newreac_gef150_k4a0_k2b'
parameters=[0.0394, 0.0788, 0.1575, 0.315]
filenames = ['dynN80k2b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[15,15,15,10]

subfolder = 'basegory_newreac_gef400_k4a0_k2b'
parameters=[1,3]
filenames = ['dynN80k2b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5]

subfolder = 'newreac_controls/gef50_k1b'
parameters=[3,30,100,300]
filenames = ['dynN80k1b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5]

subfolder = 'newreac_controls/gef100_k1b'
parameters=[3,30,100,300]
filenames = ['dynN80k1b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5]

subfolder = 'newreac_controls/gef300_k1b'
parameters=[1,3,30,100,300]
filenames = ['dynN80k1b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[15,15,5,5,5]

subfolder = 'newreac_controls/gef150_k4a_0_k1b'
parameters=[3,30,100,300]
filenames = ['dynN80k1b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5]

subfolder = 'newreac_controls/gef250_k4a_0_k1b'
parameters=[3,30,100,300]
filenames = ['dynN80k1b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5]

subfolder = 'newreac_controls/gef400_k4a_0_k1b'
parameters=[1,3,30,100,300]
filenames = ['dynN80k1b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5]

subfolder = 'newreac_controls/gef50_k5b'
parameters=[0.195,0.65,1.95,19.5,65]
filenames = ['dynN80k5b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5]

subfolder = 'newreac_controls/gef100_k5b'
parameters=[0.195,0.65,1.95,19.5,65]
filenames = ['dynN80k5b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5]

subfolder = 'newreac_controls/gef300_k5b'
parameters=[0.195,0.65,1.95,19.5,65]
filenames = ['dynN80k5b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[10,10,10,10,10]



subfolder = 'newreac_controls/gef50_k4b'
parameters=[30,100,200,500]
filenames = ['dynN80k4b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,15]

subfolder = 'newreac_controls/gef100_k4b'
parameters=[30,100,200,500]
filenames = ['dynN80k4b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5]

subfolder = 'newreac_controls/gef300_k4b'
parameters=[3,30,100,200,500]
filenames = ['dynN80k4b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[15,15,5,5,5]

subfolder = 'newreac_controls/gef150_k4a_0_k5a'
parameters=[0.04,0.12,0.4,1.2]
filenames = ['dynN80k5a'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5]

subfolder = 'newreac_controls/gef250_k4a_0_k5a'
parameters=[0.04,0.12,0.4,1.2,12]
filenames = ['dynN80k5a'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5]

subfolder = 'newreac_controls/gef400_k4a_0_k5a'
parameters=[0.04,0.12,0.4,1.2,12]
filenames = ['dynN80k5a'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5]

subfolder = 'newreac_controls/gef250_k4a_0_k5b'
parameters=[0.65,1.95]
filenames = ['dynN80k5b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5]

subfolder = 'newreac_controls/gef150_k4a_0_k5b'
parameters=[1.95,19.5]
filenames = ['dynN80k5b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5]

subfolder = 'newreac_controls/gef400_k4a_0_k5b'
parameters=[0.195,0.65,1.95]
filenames = ['dynN80k5b'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5]

#subfolder = 'basegory_newreac_gef400_k4a0_k2b'
#parameters=[1,3,10,30]
#filenames = ['dynN80k2b'+str(parameters[i])+'sim' for i in range(len(parameters))]
#skipframes=[10,5,5,5]

subfolder = 'basegory_newreac_gef100_k3_0p07_k4a'
parameters=[0.0001,0.0003,0.001,0.002,0.003,0.005,0.01,0.02,0.03,0.05,0.1,1,10]
filenames = ['dynN80k4a'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5,5,5,5,5,5,5,10,10]

subfolder = 'basegory_newreac_gef'
parameters=[15,25,50,100,200,300,500]
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,10,15,25]

subfolder = 'basegory_newreac_k4a_0_gef'
#gef 50 didnt polarize
parameters=[100,150, 200, 300,500]
filenames = ['dynN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]
skipframes=[5,5,5,5,5,5]

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
Nspecies=2
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
            
            #Get totalCdc42T = Cdc42T+GEF42 
            TS_tot42T = TS_cdc42T_gef42[0::2,:] + TS_cdc42T_gef42[1::2,:]
            TS_gef42T = TS_cdc42T_gef42[1::2,:]
            
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
            
            #GEF42
            gef42tvst = np.sum(TS_gef42T[:,:],1)
            #sum over time            
            gef42t+=np.sum(gef42tvst)            
            gef42t2+=np.sum(gef42tvst**2)
            
            ntot42t+=Ntimes            
        else:
            print "file " + name + " not found"
           
    stepshistograms.append(stepsizes)
    meantot42t.append(1.0*tot42t/ntot42t)
    meantot42t2.append(1.0*tot42t2/ntot42t)
    meangef42t.append(1.0*gef42t/ntot42t)
    meangef42t2.append(1.0*gef42t2/ntot42t)



data = {'ssd':sum_squared_disps, 'nsq':n_squared_disps, 'params':parameters,
        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
        'mgef42t': meangef42t, 'mgef42t2': meangef42t2,
        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
        'stphist':stepshistograms, 'mxstp':maxstep}

with open(folder+subfolder+'.pkl', 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open('filename.pickle', 'rb') as handle:
 #   b = pickle.load(handle)

#erase folder 
#shutil.rmtree(folder+subfolder)


