
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io
import os
import os.path
from functions_patch import get_centroid
plt.ion()


matplotlib.rcParams.update({'font.size': 12})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1.5


filenames = ['./basegory_newreac0_gef/dynN80gef700']

filenames=['./basegory_k2b_0p0625x_k5a_10x_gef_long/dynN80gef15'] 

filenames=['./basegory_k7_50_gef_long/dynN80gef15'] 

filenames=['./basegory_newreac_gef/dynN80gef15'] 


Nsims=1
sims=[0,1,2,3]
Nspecies=1
skipframes=0
centerframe=3; #use frame 4 as reference to center

for filename in filenames:
    TS_tot42T=0 #empty cdc42t distribution to avoid reusing when a .dat file is not present
    #for sim in range(Nsims):
    for sim in sims:
        #print(sim)
        name = filename + 'sim' +str(sim+1) + '.dat'
        #name =folder + filename + 'sim' + str(sim+1) + '.dat'
         
        if os.path.isfile(name):
            #print(name)
            TS_tot42T = np.genfromtxt(name,skip_header=skipframes*Nspecies)
        else:
            print "file " + name + " not found"
        
        if Nspecies != 1: 
            #get just the first species
            print "more than 1 species"
            #TS_tot42T=TS_tot42T[0::Nspecies,:]
            TS_tot42T = TS_tot42T[0::2,:] + TS_tot42T[1::2,:]
        
        Nrows=np.shape(TS_tot42T)[0]    
        Ntimes=Nrows           
        N=int(np.sqrt(np.shape(TS_tot42T)[1]))
       
        #CENTER
        #get coordinates of max value in the reference frame
        mat_first = np.reshape(TS_tot42T[centerframe,:],(N,N))
        im,jm = np.unravel_index(mat_first.argmax(), mat_first.shape)
        
        TS_centered=np.zeros((Ntimes,N*N))
        for i in range(Ntimes):
            #set as matrix
            mat=np.reshape(TS_tot42T[i,:],(N,N))
            #center
            centered=mat[np.ix_((np.arange(N) + im - int(N/2)) % N , (np.arange(N) + jm - int(N/2)) % N)]
            #store as vector
            TS_centered[i,:]=centered.flat
          
        #Save to file
        np.savetxt(filename+'centered' + str(sim+1) + '.dat', TS_centered, delimiter=' ',fmt='%d')
        
        #Compute centroid trajectory
        
        TS_centroid=[]
        for i in range(Ntimes):
            if np.sum(TS_centered[i,:]) > 0:
                mat = np.reshape(TS_centered[i,:],(N,N))
                [xcent, ycent] = get_centroid(mat)
                TS_centroid.append([xcent, ycent])
            else:
                TS_centroid.append([np.nan, np.nan])
                
        TS_centroid=np.asarray(TS_centroid)    
        np.savetxt(filename+'centroid' + str(sim+1) + '.dat', TS_centroid, delimiter=' ')
        

    #plt.hist(TS_centered[30,:],bins=np.arange(0,20,0.5))
    #plt.ylim([0,100])

