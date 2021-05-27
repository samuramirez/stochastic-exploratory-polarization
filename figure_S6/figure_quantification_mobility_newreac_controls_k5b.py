
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



folder='D:/dynamic_polarity_data/newreac_controls/'


subfolder = 'gef50_k5b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data = pickle.load(handle)
    
subfolder = 'gef100_k5b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data2 = pickle.load(handle)

subfolder = 'gef300_k5b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data3 = pickle.load(handle)

subfolder = 'gef150_k4a_0_k5b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data4 = pickle.load(handle)

subfolder = 'gef250_k4a_0_k5b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data5 = pickle.load(handle)

subfolder = 'gef400_k4a_0_k5b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data6 = pickle.load(handle)


    
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
nlinear=12 #?
n0=4   
for i in range(len(data['params'])):
    msd=data['ssd'][i,:]/data['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data['smpl']*(np.arange(len(msd)))/60.    
    #plt.plot(intervals,msd,label='D='+'{:.4f}'.format(Dpatch) + r' ($\mu m^2/min$)', color=colors[i])
    p,cov = np.polyfit(np.log(intervals[n0:nlinear+n0]),np.log(msd[n0:nlinear+n0]),1,cov='True')
    Dp.append(np.exp(p[1])/4.0)
    a.append(p[0])
    stdDp.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='black',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    #plt.ylim([0,1])
    plt.xlim([0,60])
    plt.yscale('log')
    plt.xscale('log')
    #plt.legend(loc=2,fontsize=12)
    fig.set_size_inches(2.5,1.77)   
#plt.title('Max step considered='+str(maxstep))
plt.tight_layout()


#MSD 
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp2=[]
a2=[]
stdDp2=[]
stda2=[]
nlinear=12 #?
n0=4   
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

#MSD 
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp3=[]
a3=[]
stdDp3=[]
stda3=[]
nlinear=12 #?
n0=4   
for i in range(len(data3['params'])):
    msd=data3['ssd'][i,:]/data3['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data3['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[n0:nlinear+n0]),np.log(msd[n0:nlinear+n0]),1,cov='True')
    Dp3.append(np.exp(p[1])/4.0)
    a3.append(p[0])
    stdDp3.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda3.append(np.sqrt(np.diag(cov)[0]))        
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

#MSD 
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp4=[]
a4=[]
stdDp4=[]
stda4=[]
n0=1
nlinear=10 #?    
for i in range(len(data4['params'])):
    msd=data4['ssd'][i,:]/data4['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data4['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[n0:nlinear+n0]),np.log(msd[n0:nlinear+n0]),1,cov='True')
    Dp4.append(np.exp(p[1])/4.0)
    a4.append(p[0])
    stdDp4.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda4.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='black',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,30])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()

#MSD 
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp5=[]
a5=[]
stdDp5=[]
stda5=[]
n0=1
nlinear=10 #?    
for i in range(len(data5['params'])):
    msd=data5['ssd'][i,:]/data5['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data5['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[n0:nlinear+n0]),np.log(msd[n0:nlinear+n0]),1,cov='True')
    Dp5.append(np.exp(p[1])/4.0)
    a5.append(p[0])
    stdDp5.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda5.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='black',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,30])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()

#MSD 
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp6=[]
a6=[]
stdDp6=[]
stda6=[]
n0=4
nlinear=12 #?    
for i in range(len(data6['params'])):
    msd=data6['ssd'][i,:]/data6['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data6['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[n0:nlinear+n0]),np.log(msd[n0:nlinear+n0]),1,cov='True')
    Dp6.append(np.exp(p[1])/4.0)
    a6.append(p[0])
    stdDp6.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda6.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='black',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,30])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()



plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots() 
ax.errorbar(data['params'],Dp,stdDp,color='black',marker="o",markersize=2,capsize=2,label=r'GEF=50',linestyle='dotted')
ax.errorbar(data2['params'],Dp2,stdDp2,color='black',marker="o",markersize=2,capsize=2,label=r'GEF=100',linestyle='--')
ax.errorbar(data3['params'],Dp3,stdDp3,color='black',marker="o",markersize=2,capsize=2,label=r'GEF=300',linestyle='-')
ax.set_xlabel(r'$k_{5b}(1/s)$')
ax.set_title(r'$k_{4a}=2\mu m^2/s$')
#ax.set_xticks([0,100,200,300,400,500])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
ax.legend(loc=2,fontsize=5)
fig.set_size_inches(2.3,2)
ax.set_ylim([5e-4,10])
#plt.axhline(y=1./60,color='black',linestyle='--')
#ax.set_xlim([0,510])
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
fig.savefig('Dpatch_newreac_controls_k5b_gef.pdf')
plt.show()


plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots() 
ax.errorbar(data4['params'],Dp4,stdDp4,color='black',marker="o",markersize=2,capsize=2,label=r'GEF=150',linestyle='dotted')
ax.errorbar(data5['params'],Dp5,stdDp5,color='black',marker="o",markersize=2,capsize=2,label=r'GEF=250',linestyle='--')
ax.errorbar(data6['params'],Dp6,stdDp6,color='black',marker="o",markersize=2,capsize=2,label=r'GEF=400',linestyle='-')
ax.set_xlabel(r'$k_{5b} (1/s)$')
ax.set_title(r'$k_{4a}=0$')
#ax.set_xticks([0,100,200,300,400,500])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
ax.legend(loc=2,fontsize=5)
fig.set_size_inches(2.3,2)
ax.set_ylim([5e-4,10])
#plt.axhline(y=1./60,color='black',linestyle='--')
#ax.set_xlim([0,510])
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
fig.savefig('Dpatch_newreac_controls_k5b_k4a_0_gef.pdf')
plt.show()

