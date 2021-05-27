
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


#folder='./newreac_gef100_k4a_k3/'
#subfolder='basegory_newreac_gef100_k3_0p07_k4a'
#parameters=[1,0.1, 0.01,0.005,0.002,0]
#parameters=[0,0.002, 0.005 ]
#filenames = ['dynN80k4a'+str(parameter)+'sim' for parameter in parameters ]


folder = './'
subfolder = 'basegory_newreac_gef100_k3_0p07_k4a'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data = pickle.load(handle)


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
nlinear=8
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
matplotlib.rcParams.update({'font.size': 8})
fig,ax = plt.subplots() 
ax.errorbar(data['params'],Dp,stdDp,color='black',marker="o",markersize=2,capsize=2)
ax.set_xlabel(r'$k_{4a}$ ($\mu m^2/s$)')
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
#plt.legend(loc=1,fontsize=9)
fig.set_size_inches(3,2.5)
#ax.set_ylim([0,0.25])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([3e-4,20])
ax.set_ylim([1e-3,2])
ax.set_xscale('log')
ax.set_yscale('log')

#ax2=ax.twinx()
# make a plot with different y-axis using second axis object
#stddevtot42t= (np.asarray(meantot42t2) - np.asarray(meantot42t)**2)**0.5
#stderrtot42t= stddevtot42t/((Ntimes-5)*Nsims)**0.5
#ax2.errorbar(parameters, meantot42t,stderrtot42t,color="red",marker="o",markersize=2,capsize=2)
#ax2.scatter(parameters,Dp,color='blue',s=5)
#ax2.set_ylabel("Total Cdc42T")
#ax2.spines['right'].set_color('red')
#ax2.tick_params(axis='y', colors='red')
#ax2.set_yticks([0,500,1000,1500,2000,2500])
#ax2.yaxis.label.set_color('red')
#ax2.set_ylim([0,2700])

plt.tight_layout()
fig.savefig('Dpatch_newreac_gef100_k4a.pdf')

#SAVE DATA
dataD = {'params':data['params'], 'Dp':Dp, 'stdDp':stdDp, 'a':a, 'stda':stda}
with open(folder+subfolder+'_dataD.pkl', 'wb') as handle:
    pickle.dump(dataD, handle, protocol=pickle.HIGHEST_PROTOCOL)