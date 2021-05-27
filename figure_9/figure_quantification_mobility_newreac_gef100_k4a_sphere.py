
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

import scipy.io
import os
import os.path
from scipy.optimize import curve_fit
from functions_patch import get_centroid
from functions_patch import sum_squared_displacements
import zipfile

plt.ion()



folder='D:/dynamic_polarity_data/fig_3D_transition_quant/'
subfolder = 'newreac_gef100_k4a'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data = pickle.load(handle)

#here data['params'] are saved as strings, below I make them numbers
data['params']=[0.001,0.002,0.005,0.01,0.02,0.05,0.1,1]


folder = 'D:dynamic_polarity_data/'
subfolder = 'basegory_newreac_gef100_k3_0p07_k4a'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data2 = pickle.load(handle)


#labels=[str(parameters[i]) for i in range(len(parameters))]
#colors = [cm.jet(i) for i in np.linspace(1, 0, len(parameters))]

#PLOTS     
matplotlib.rcParams.update({'font.size': 8})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1


plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp=[]
a=[]
stdDp=[]
stda=[]
n0=1
nlinear=6
for i in range(len(data['params'])):
    msd=data['ssd'][i,:]/data['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[n0:nlinear+n0]),np.log(msd[n0:nlinear+n0]),1,cov='True')
    Dp.append(np.exp(p[1])/4.0)
    a.append(p[0])
    stdDp.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='black',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([1e-4,20])
    plt.xlim([0.9,25])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()

plt.close('all')
fig=plt.figure(fi)
fi+=1
Dp2=[]
a2=[]
stdDp2=[]
stda2=[]
n0=1
nlinear=8
for i in range(len(data2['params'])):
    msd=data2['ssd'][i,:]/data2['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data2['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[n0:nlinear+n0]),np.log(msd[n0:nlinear+n0]),1,cov='True')
    Dp2.append(np.exp(p[1])/4.0)
    a2.append(p[0])
    stdDp2.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda2.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='black',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([1e-4,20])
    plt.xlim([0.9,25])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()



plt.close('all')
matplotlib.rcParams.update({'font.size': 8})
fig,ax = plt.subplots() 
parameters = data['params']
#ax.plot(parameters,Dp,color='black',marker="o",markersize=3)
ax.errorbar(parameters,Dp,stdDp,color='red',marker="o",markersize=2,capsize=2,label='3D particle-based')
ax.errorbar(data2['params'],Dp2,stdDp2,color='black',marker="o",markersize=2,capsize=2,label='2D Spatial Gillespie')

ax.set_xlabel(r'$k_{4a}$ ($\mu m^2/s$)')
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
plt.legend(loc=1,fontsize=6)
fig.set_size_inches(3,2.5)
#ax.set_ylim([0,0.25])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([3e-4,20])
ax.set_ylim([1e-3,2])
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
fig.savefig('Dpatch_newreac_gef100_k4a_sphere.pdf')

