
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io
import os
import os.path
import zipfile
from scipy.optimize import curve_fit
from functions_patch import get_centroid
from functions_patch import sum_squared_displacements
import pandas as pd
import shutil
import pickle


folder='D:/dynamic_polarity_data/newreac_mm_gef/'

parameters=[0.0001,0.001,0.01,0.1,1]

filenames = ['dynN80k3'+str(parameter)+'sim' for parameter in parameters ]
subfolders=['basegory_newreac_gef100_k2_0_k4a_0_k3','basegory_newreac_gef200_k2_0_k4a_0_k3','basegory_newreac_gef300_k2_0_k4a_0_k3','basegory_newreac_gef400_k2_0_k4a_0_k3']

filenames = ['dynN80k4a'+str(parameter)+'sim' for parameter in parameters ]
subfolders=['basegory_newreac_gef100_k2_0_k3_0_k4a','basegory_newreac_gef200_k2_0_k3_0_k4a','basegory_newreac_gef300_k2_0_k3_0_k4a','basegory_newreac_gef400_k2_0_k3_0_k4a']

filenames = ['dynN80k2a'+str(parameter)+'sim' for parameter in parameters ]
subfolders=['basegory_newreac_gef100_k3_0_k4a_0_k2a','basegory_newreac_gef200_k3_0_k4a_0_k2a','basegory_newreac_gef300_k3_0_k4a_0_k2a','basegory_newreac_gef400_k3_0_k4a_0_k2a']


for subfolder in subfolders:    
    #extract files
    if not os.path.isdir(folder+subfolder):
        with zipfile.ZipFile(folder+subfolder + '.zip',"r") as zip_ref:
            zip_ref.extractall(folder)
    
    N1=80
    Nsims=50
    T=1*3600
    sampling= 1 #1min intervals for centroid calculation
    sampling_sim = 60 #s sampling time in simulations
    skipframes=5 #first point analyzed is 1min
    interval=T/sampling_sim #nframes
    nframes=T/sampling_sim - skipframes
    samples=int(nframes/sampling)
    Nspecies=2
    L=8.0
    centerwrt=1 #center with respect to frame

       
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
        
        for sim in range(Nsims):
            print(sim )
            name =folder + subfolder + '/' + filename + str(sim+1) + '.dat'
          
            if os.path.isfile(name):
                #print(name)
                TS_cdc42T_gef42 = np.genfromtxt(name,skip_header=skipframes*Nspecies)       
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
            
        meantot42t.append(1.0*tot42t/ntot42t)
        meantot42t2.append(1.0*tot42t2/ntot42t)
        meangef42t.append(1.0*gef42t/ntot42t)
        meangef42t2.append(1.0*gef42t2/ntot42t)
       
    Dp=[]
    stdDp=[]
    a=[]
    stda=[]
    nlinear=7 #?   
    n0=1
    for i in range(len(parameters)):
        msd=sum_squared_disps[i,:]/n_squared_disps[i,:]
        msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
        intervals=sampling*sampling_sim*(np.arange(len(msd)))/60.
        p,cov = np.polyfit(np.log(intervals[n0:nlinear+n0]),np.log(msd[n0:nlinear+n0]),1,cov='True')
        Dp.append(np.exp(p[1])/4.0)
        a.append(p[0])
        stdDp.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
        stda.append(np.sqrt(np.diag(cov)[0]))
        
       
    #Save the data
    data = pd.DataFrame([Dp,stdDp,a,stda])
    data.columns=parameters
    data.index=['Dpatch','stdDpatch','a','stda']
    data.to_pickle(folder + subfolder+'.pkl')
    
    alldata = {'ssd':sum_squared_disps, 'nsq':n_squared_disps, 'params':parameters,
        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
        'mgef42t': meangef42t, 'mgef42t2': meangef42t2,
        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim}

    with open(folder+subfolder+'_all.pkl', 'wb') as handle:
        pickle.dump(alldata, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    #erase folder 
    shutil.rmtree(folder+subfolder)
        
    #to read 
    #data = pd.read_pickle(folder + subfolder+'.pkl')
    

