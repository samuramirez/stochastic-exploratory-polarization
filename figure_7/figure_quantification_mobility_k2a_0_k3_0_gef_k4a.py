
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io
import os
import os.path
from scipy.optimize import curve_fit
import pandas as pd

plt.ion()


matplotlib.rcParams.update({'font.size': 8})
#matplotlib.rcParams.update({'linewidth': 4})

matplotlib.rcParams['lines.linewidth'] = 1.5

'''
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}

matplotlib.rc('font', **font)
'''

allpoints=[]

folder='./newreac_mm_gef/'
figurename='figure_mobility_Dpatch_gef_vs_k4a.pdf'

metric='Dpatch'  
factsize=300

subfolder='basegory_newreac_gef100_k2_0_k3_0_k4a'
data = pd.read_pickle(folder + subfolder+'.pkl')
gef=100;
points=[ [gef, k4a, data.loc[metric][k4a]]  for k4a in data.columns ] 
allpoints=allpoints + points

subfolder='basegory_newreac_gef200_k2_0_k3_0_k4a'
data = pd.read_pickle(folder + subfolder+'.pkl')
gef=200;
points=[ [gef, k4a, data.loc[metric][k4a]]  for k4a in data.columns ] 
allpoints=allpoints + points

subfolder='basegory_newreac_gef300_k2_0_k3_0_k4a'
data = pd.read_pickle(folder + subfolder+'.pkl')
gef=300;
points=[ [gef, k4a, data.loc[metric][k4a]]  for k4a in data.columns ] 
allpoints=allpoints + points

subfolder='basegory_newreac_gef400_k2_0_k3_0_k4a'
data = pd.read_pickle(folder + subfolder+'.pkl')
gef=400;
points=[ [gef, k4a, data.loc[metric][k4a]]  for k4a in data.columns ] 
allpoints=allpoints + points


plt.close('all')

fig=plt.figure(1)
for i in range(len(allpoints)):
    plt.scatter(allpoints[i][1],allpoints[i][0], s=factsize*allpoints[i][2], c='black' )




plt.xscale('log')
plt.yscale('linear')

#plt.xlim([-0.0005 ,1.1e-2])
plt.xlim([5e-4 ,2])
plt.ylim([0,450])
plt.ylabel(r'GEF')
plt.xlabel(r'$k_{4a}$ ($\mu m^2/s\min$)')
#plt.title('Mean step size\n(1 min intervals)',fontsize=8)
plt.title(r'$D_{patch}$ ($\mu m^2/min$)',fontsize=8)

#for legend
#plt.scatter(10,300, s=factsize*0.2, c='black' )
#plt.scatter(10,250, s=factsize*0.05, c='black' )
#plt.scatter(10,200, s=factsize*0.005, c='black' )
#plt.xlim([5e-4 ,20])

fig.set_size_inches(3, 2.5)
plt.tight_layout()
fig.savefig(figurename)


'''

#PLOTS
plt.close('all')
fi=1  
#MSD
fig=plt.figure(fi)
fi+=1
Dp=[]
nlinear=10 #?    
for i in range(len(parameters)):
    msd=sum_squared_disps[i,:]/n_squared_disps[i,:]
    msd=np.concatenate(([0],msd))
    intervals=sampling*sampling_sim*(np.arange(len(msd)))/60.
    Dpatch = np.polyfit(intervals[:nlinear],msd[:nlinear],1)[0]/4
    Dp.append(Dpatch)
    plt.plot(intervals,msd,label='D='+'{:.4f}'.format(Dpatch) + r' ($\mu m^2/min$)', color=colors[i])
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.legend(loc=1,fontsize=9)
    fig.set_size_inches(3,2.5)
    plt.ylim([0,15])
    plt.xlim([0,30])
#plt.title('Max step considered='+str(maxstep))
plt.tight_layout()

#D vs parameter
fig=plt.figure(fi)
fi+=1
plt.plot(parameters,Dp)
plt.xlabel(r'$k_{4a}$ ($\mu m^2/s$)')
plt.ylabel(r'Deff ($\mu m^2/min$)')
#plt.legend(loc=1,fontsize=9)
fig.set_size_inches(3,2.5)
#plt.ylim([0,10])
#plt.xlim([0,30])
#plt.title('Max step considered='+str(maxstep))
plt.tight_layout()

#STEPS HISTOGRAM
meansteps=[]
stdsteps=[]
nsteps=[]
for i in range(len(parameters)):
    fig=plt.figure(fi)
    fi+=1
    #bins=np.arange(0.0, 40, 1)
    #keep non nan values
    steps = stepshistograms[i][~np.isnan(stepshistograms[i])]
    plt.hist(steps,bins=40,color='blue',density=True,alpha=1,label= labels[i]  + ' n='+str(len(steps)) ,histtype='step',linewidth=2)
    #plt.hist(stepsizesgt0,bins=40,color='blue',density=True,alpha=1,label='N ='+str(N)+' Dc=5 n='+str(len(stepsizesgt0)) ,histtype='step',linewidth=2)
    meanstep=np.nanmean(stepshistograms[i])
    meansteps.append(meanstep)
    stdstep=np.nanstd(stepshistograms[i])
    stdsteps.append(stdstep)
    nsteps.append(np.sum(~np.isnan(stepshistograms[i])))
    plt.axvline(x=meanstep,color='blue')
    plt.xlim([0 , 5])
    plt.xlabel(r'Step size ($\mu$m)')
    plt.ylabel('Density')
    fig.set_size_inches(5, 4)
    plt.legend()
    plt.tight_layout()

#meanstep vs parameter
fig=plt.figure(fi)
fi+=1
plt.plot(parameters,meansteps)
plt.xlabel(r'$k_{4a}$ ($\mu m^2/s$)')
plt.ylabel(r'Mean step ($\mu m$)')
#plt.legend(loc=1,fontsize=9)
fig.set_size_inches(3,2.5)
#plt.ylim([0,10])
#plt.xlim([0,30])
#plt.title('Max step considered='+str(maxstep))
plt.tight_layout()




'''