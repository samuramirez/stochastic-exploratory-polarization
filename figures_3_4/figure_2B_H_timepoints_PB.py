import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io
import os.path
from matplotlib import cm
import pandas as pd
plt.ion()
matplotlib.rcParams.update({'font.size': 18})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1.5

Nsims=30;
GEFconc=["200","300","400","500","600","700"];

#GEFconc=["500","600","700"];
#timeanalysis=[1:Nf];%verify this with H function file
#for i =1:length(GEFconc)

#Hallsims(l,:,:)=H_vals; %nsims x 50Hvalues * 30timepoints 

#verify this with matlab computation of H
max_r  = 5
dr=0.1
#in matlab r_vals = 0:dr:max_r;
r_vals = np.arange(0,max_r+dr,dr)
times=np.arange(10,310,10)
L=8
protr=0.02
#meth=['ka','kh','ksrhet'];
meth=['kh','ksrhet'];
meth=['ksrhet'];

Nmeth=["40","80","160"];
Nmeth=["120"];

methodname=[r'$k_h$',r'$k_c$']
tpoints=np.asarray([10,60, 120, 300])/10 -1;
sim=0
gefconc='700'
colors = [cm.jet(i) for i in np.linspace(1, 0, len(GEFconc))]
#colNs= [cm.jet(i) for i in np.linspace(1, 0, len(Nmeth))] 
cols= [ cm.jet(0.2) ,'red']
#markers=['x','^']
linestyles=['dotted','dashed','solid']

fi=1
fig=plt.figure(fi)
fi+=1
fig, axes = plt.subplots(nrows=1, ncols=len(tpoints), figsize=(12, 3))
for idxtpoints,tpoint in enumerate(tpoints):
    
    #PARTICLE BASED
    #Maximum of H function in particle-based simulations is consistently r=1.1
    file =  filename='./HvaluesPB/Hvaluesgef' + gefconc + 'PB.mat';
    #file =  filename='Hvaluesgef' + gefconc + 'PB.mat';

    if os.path.isfile(file):
        matdata = scipy.io.loadmat(file)
        
        Hallsims=matdata['Hallsims']
        #look at individual simulations
        
        #plt.plot(r_vals,Hallsims[sim,:,tpoint])
        #average the last timepoints to account for noisy clustering
        axes[idxtpoints].plot(r_vals,np.nanmean(Hallsims[sim,:,tpoint:tpoint+2],1),color='black')
        
        #look at average of individual simulations    
        #plt.plot(r_vals,np.nanmean(Hallsims[:,:,25:30],(0,2)),color=colors[idxgef],label=GEFconc[idxgef])
        #plt.plot(r_vals,np.nanmean( np.nanmean(Hallsims[:,:,25:30],2) ,0))    
    
        #plt.plot(r_vals,np.mean(Hallsims[:,:,25],0))
        #plt.title(metstitle[imet])
        axes[idxtpoints].set_xlabel(r'r ($\mu m$)')
        #plt.ylabel('H(r)')
        axes[idxtpoints].set_ylim(-0.5,3.5)
        axes[idxtpoints].set_xlim(-0.5,4.5)
        axes[idxtpoints].set_xticks([0,1,2,3,4])
        axes[idxtpoints].set_yticks([0,1.5,3])

        #fig.set_size_inches(5, 4)
        #fig.set_size_inches(3, 2.5)
        #plt.legend(fontsize=8)
    else:
        print("file " + file + " not found")
plt.tight_layout()
fig.savefig('H_PB_timepoints.pdf')

'''   
    #SPATIAL GILLESPIE    
    file =  filename='./Hvalues/Hvalues'+'gef'+ gefconc +meth[0]+'N'+Nmeth[0]+'.mat';

    if os.path.isfile(file):
        matdata = scipy.io.loadmat(file)        
        Hallsims=matdata['Hallsims']      
        #look at average of individual simulations    
        #plt.plot(r_vals,np.nanmean(Hallsims[:,:,25:30],(0,2)),color=colors[idxgef],linestyle='--')
        plt.xlabel(r'separation ($\mu m$)')
        plt.ylabel('H')
        plt.ylim(-0.5,3.5)
        #plt.xlim(0,1.3)
        fig.set_size_inches(5, 4)
        #plt.legend(fontsize=8) 
        plt.tight_layout()
    else:
        print("file " + file + " not found")
'''
