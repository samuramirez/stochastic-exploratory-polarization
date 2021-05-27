
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



folder = './'
subfolder = 'basegory_newreac_gef'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data = pickle.load(handle)
    
  
#data = {'ssd':sum_squared_disps, 'nsq':n_squared_disps, 'params':parameters,
#        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
#        'mgef42t': meangef42t, 'mgef42t2': meangef42t2,
#        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
#        'stphist':stepshistograms}

matplotlib.rcParams.update({'font.size': 7})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1


#PLOTS
 #colors=['black','magenta','green','cyan']
#colors=['black','grey','green','lightgrey']
#linestyles=['-','--','dotted','dotted']
#colors=[cm.Greys(100),cm.Greys(50),cm.Greys(50),cm.Greys(10)]
    
labels= [str(data['params'][i]) for i in range(len(data['params']))]
#colors=['blue','red','green','black']
#colors = [ cm.jet(0.2), cm.jet(0.5) , cm.jet(0.7) ,cm.jet(0.9) ]
colors = [cm.hot(i) for i in np.linspace(0.5, 0, len(data['params']))]

matplotlib.rcParams.update({'font.size': 7})


#MSD 

plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp=[]
a=[]
stdDp=[]
stda=[]
for i in range(len(data['params'])):
    if data['params'][i] >= 200  :
        n0=3
        nlinear=12
    else:
        n0=1
        nlinear=8
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



#THIS IS FIGURE 5E.  Dpatch vs GEF for k4a=2 and k4a=0 LOG
plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots() 
#ax.plot(parameters,Dp,color='black',marker="o",markersize=3)
ax.errorbar(data['params'],Dp,stdDp,color='black',marker="o",markersize=2,capsize=2)
#ax.errorbar(data2['params'],Dp2,stdDp2,color='red',marker="o",markersize=2,capsize=2,label=r'$k_{4a}=0$, $k_{7}=0.5\mu m^2/s$ ')
#ax.errorbar(data3['params'],Dp3,stdDp3,color='blue',marker="o",markersize=2,capsize=2,label=r'$k_{4a}=0$, $k_{7}=0.5\mu m^2/s$ ')
#ax.errorbar(data4['params'],Dp4,stdDp4,color='magenta',marker="o",markersize=2,capsize=2,label=r'$k_{4a}=0$, $k_{7}=0.5\mu m^2/s$ ')

#ax2=ax.twinx()
# make a plot with different y-axis using second axis object
#stddevgef42t= (np.asarray(data['mgef42t2']) - np.asarray(data['mgef42t'])**2)**0.5
#ax2.errorbar(data['params'],data['mgef42t'],stddevgef42t,color='orange',marker="o",markersize=2,capsize=2,linestyle='-')
#ax2.set_ylabel("Cdc42T-GEF")
#ax2.spines['right'].set_color('orange')
#ax2.tick_params(axis='y', colors='orange')
#ax2.set_yticks([0,500,1000,1500,2000,2500])
#ax2.yaxis.label.set_color('orange')
#ax2.set_ylim([0,500])

ax.set_xlabel(r'GEF')
ax.set_xticks([0,100,200,300,400,500])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
#ax.legend(loc=1,fontsize=7)
fig.set_size_inches(2.3,1.77)
#fig.set_size_inches(2.5,1.77) # after adding second axis

ax.set_ylim([5e-4,1])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([0,510])
#ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
#fig.savefig('Dpatch_newreac_gef.pdf')
plt.show()

#SAVE DATA
dataD = {'params':data['params'], 'Dp':Dp, 'stdDp':stdDp, 'a':a, 'stda':stda}
with open(folder+subfolder+'_dataD.pkl', 'wb') as handle:
    pickle.dump(dataD, handle, protocol=pickle.HIGHEST_PROTOCOL)


