
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io
import os
import os.path
from scipy.optimize import curve_fit
import pandas as pd
import pickle

plt.ion()


matplotlib.rcParams.update({'font.size': 8})
#matplotlib.rcParams.update({'linewidth': 4})

matplotlib.rcParams['lines.linewidth'] = 1.0

#alldata = {'ssd':sum_squared_disps, 'nsq':n_squared_disps, 'params':parameters,
#        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
#        'mgef42t': meangef42t, 'mgef42t2': meangef42t2,
#        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
#        'stphist':stepshistograms}

'''
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}

matplotlib.rc('font', **font)
'''


folder='./newreac_mm_gef/'


figurename='figure_mobility_Dpatch_gef_vs_k3.pdf'
subfsufix='_k2_0_k4a_0_k3'
xlabel=r'$k_{3a}$ ($\mu m^2/min$)'

#figurename='figure_mobility_Dpatch_gef_vs_k2a_legend.pdf'
#subfsufix='_k3_0_k4a_0_k2a'
#xlabel=r'$k_{2a}$ ($\mu m^2/min$)'

#figurename='figure_mobility_Dpatch_gef_vs_k4a.pdf'
#subfsufix='_k2_0_k3_0_k4a'
#xlabel=r'$k_{4a}$ ($\mu m^2/min$)'

metric='Dpatch'  
factsize=300

allpoints=[]
tot42tvsgef=[]
gef42tvsgef=[] #indices: gef, mean or mean^2, k 

gefs=[100,200,300,400]
for gef in gefs:
    subfolder='basegory_newreac_gef'+str(gef)+subfsufix
    #read data
    data = pd.read_pickle(folder + subfolder+'.pkl')
    with open(folder+subfolder+'_all.pkl', 'rb') as handle:
        alldata = pickle.load(handle)    
    
    points=[ [gef, k, data[k][metric]]  for k in data.columns ] 
    allpoints=allpoints + points
    tot42tvsgef.append([alldata['mtot42t'], alldata['mtot42t2'] ])
    gef42tvsgef.append([alldata['mgef42t'], alldata['mgef42t2'] ])
    



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
plt.xlabel(xlabel)
#plt.title('Mean step size\n(1 min intervals)',fontsize=8)
plt.title(r'$D_{patch}$ ($\mu m^2/min$)',fontsize=8)

#for legend
#plt.scatter(10,300, s=factsize*0.2, c='black' )
#plt.scatter(10,250, s=factsize*0.05, c='black' )
#plt.scatter(10,200, s=factsize*0.005, c='black' )
#plt.xlim([5e-4 ,20])

fig.set_size_inches(3, 2.5)
plt.tight_layout()
#fig.savefig(figurename)



gef42tvsgef=np.asarray(gef42tvsgef)
#first index:gef, second index: species, last index k3

#mean totgef42t, k3a=1e-2
#massvsgef[:,0,1]


#plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots()

 #indices: gef, mean or mean^2, k 
ax.errorbar(gefs,gef42tvsgef[:,0,2],(gef42tvsgef[:,1,2] - gef42tvsgef[:,0,2]**2)**0.5,label=r'$k_3 = 10^{-2}$',color="black",linestyle='-',marker="o",markersize=2,capsize=2)
ax.errorbar(gefs,gef42tvsgef[:,0,3],(gef42tvsgef[:,1,3] - gef42tvsgef[:,0,3]**2)**0.5,label=r'$k_3 = 10^{-1}$',color="red",linestyle='-',marker="o",markersize=2,capsize=2)
ax.set_ylabel(r'Cdc42T-GEF')


ax.set_xlabel(r'GEF')

ax.legend(loc=4,fontsize=7)
fig.set_size_inches(2.3,1.77)

plt.tight_layout()
#fig.savefig('gef42t_k4a_0_k2a_0_k3_GEF.pdf')
plt.show()












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