
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



folder='./newreac_controls/'


subfolder = 'basegory_newreac_gef300_k2b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data = pickle.load(handle)
 
subfolder = 'basegory_newreac_gef200_k2b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data2 = pickle.load(handle)    

subfolder = 'basegory_newreac_gef100_k2b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data3 = pickle.load(handle)  

subfolder = 'basegory_newreac_gef50_k2b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data4 = pickle.load(handle)  


#didn't polarize with the current parameters
#subfolder = 'basegory_newreac_gef150_k4a0_k2b'
#with open(folder+subfolder+'.pkl', 'rb') as handle:
#    data5 = pickle.load(handle)  

subfolder = 'basegory_newreac_gef200_k4a0_k2b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data6 = pickle.load(handle)  

subfolder = 'basegory_newreac_gef400_k4a0_k2b'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data7 = pickle.load(handle)  



    
#data = {'ssd':sum_squared_disps, 'nsq':n_squared_disps, 'params':parameters,
#        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
#        'mgef42t': meangef42t, 'mgef42t2': meangef42t2,
#        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
#        'stphist':stepshistograms, 'mxstp':maxstep}

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

#MSD gef = 300
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp=[]
a=[]
stdDp=[]
stda=[]
nlinear=20 #?  
n0=4  
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
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='red',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,30])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()


#MSD gef = 200
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp2=[]
a2=[]
stdDp2=[]
stda2=[]
nlinear=10 #?    
for i in range(len(data2['params'])):
    msd=data2['ssd'][i,:]/data2['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data2['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[1:nlinear+1]),np.log(msd[1:nlinear+1]),1,cov='True')
    Dp2.append(np.exp(p[1])/4.0)
    a2.append(p[0])
    stdDp2.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda2.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='red',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,15])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()

#MSD gef = 100
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp3=[]
a3=[]
stdDp3=[]
stda3=[]
nlinear=10 #?    
for i in range(len(data3['params'][:-1])):
    msd=data3['ssd'][i,:]/data3['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data3['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[1:nlinear+1]),np.log(msd[1:nlinear+1]),1,cov='True')
    Dp3.append(np.exp(p[1])/4.0)
    a3.append(p[0])
    stdDp3.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda3.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='red',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,15])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()

#MSD gef = 50
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp4=[]
a4=[]
stdDp4=[]
stda4=[]
nlinear=10 #?    
for i in range(len(data4['params'])):
    msd=data4['ssd'][i,:]/data4['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data4['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[1:nlinear+1]),np.log(msd[1:nlinear+1]),1,cov='True')
    Dp4.append(np.exp(p[1])/4.0)
    a4.append(p[0])
    stdDp4.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda4.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='red',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,20])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()


plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp5=[]
a5=[]
stdDp5=[]
stda5=[]
nlinear=10 #?    
for i in range(len(data5['params'])):
    msd=data5['ssd'][i,:]/data5['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    #msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data5['smpl']*(np.arange(len(msd)))/60.    
    #p,cov = np.polyfit(np.log(intervals[1:nlinear+1]),np.log(msd[1:nlinear+1]),1,cov='True')
    #Dp5.append(np.exp(p[1])/4.0)
    #a5.append(p[0])
    #stdDp5.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    #stda5.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    #plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='red',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,30])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()

plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp6=[]
a6=[]
stdDp6=[]
stda6=[]
nlinear=15 #? 
n0=4   
#for i in range(len(data6['params'])):
for i in [2,3]:
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
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='red',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,30])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()


plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp7=[]
a7=[]
stdDp7=[]
stda7=[]
nlinear=10 #?    
#for i in range(len(data7['params'])):
for i in [0]:
    msd=data7['ssd'][i,:]/data7['nsq'][i,:]
    #msd=np.concatenate(([0],msd))
    msd=msd - msd[0] #subtract instantaneous intrinsic fluctuations of the patch
    intervals=data7['smpl']*(np.arange(len(msd)))/60.    
    p,cov = np.polyfit(np.log(intervals[1:nlinear+1]),np.log(msd[1:nlinear+1]),1,cov='True')
    Dp7.append(np.exp(p[1])/4.0)
    a7.append(p[0])
    stdDp7.append( np.exp(p[1])*np.sqrt(np.diag(cov)[1]) )
    stda7.append(np.sqrt(np.diag(cov)[0]))        
    plt.plot(intervals,msd,color='red',linewidth=1.5,marker='o')
    plt.plot(intervals,np.exp(p[1])*intervals**p[0],color='red',linewidth=1,linestyle='--')  
    plt.xlabel('Interval (min)')
    plt.ylabel(r'MSD ($\mu m^2$)')
    plt.ylim([4e-3,20])
    plt.xlim([0.9,20])
    plt.yscale('log')
    plt.xscale('log')
    fig.set_size_inches(2.5,1.77)   
plt.tight_layout()

'''
#THIS IS FIGURE 5E.  Dpatch vs GEF for k4a=2 and k4a=0
plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots() 
#ax.plot(parameters,Dp,color='black',marker="o",markersize=3)
ax.errorbar(data['params'],Dp,stdDp,color='black',marker="o",markersize=2,capsize=2,label=r'$k_{4a}=2$, $k_{7}=0.5\mu m^2/s$ ')

ax.set_xlabel(r'GEF')
ax.set_xticks([0,50,100,150,200,250,300])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
#ax.legend(loc=1,fontsize=7)
fig.set_size_inches(2.3,1.77)
#ax.set_ylim([1e-4,0.25])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([0,400])
#ax.set_xscale('log')
ax.set_yscale('linear')
plt.tight_layout()
#fig.savefig('Dpatch_newreac_k4a_2_k4a_0_gef.pdf')
plt.show()

'''
#
plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots() 
#ax.plot(parameters,Dp,color='black',marker="o",markersize=3)

ax.errorbar(data4['params'],Dp4,stdDp4,color='black',marker="o",markersize=2,capsize=2,label='GEF=50',linestyle='dotted')
ax.errorbar(data3['params'][:-1],Dp3,stdDp3,color='black',marker="o",markersize=2,capsize=2,label='GEF=100',linestyle='--')
#ax.errorbar(data2['params'],Dp2,stdDp2,color='red',marker="o",markersize=2,capsize=2,label='GEF=200')
ax.errorbar(data['params'],Dp,stdDp,color='black',marker="o",markersize=2,capsize=2,label='GEF=300',linestyle='-')
ax.set_title(r'$k_{4a}=2\mu m^2/s$')
ax.set_xlabel(r'$k_{2b}$ $(s^{-1})$')
ax.set_xticks([1,3,10,30])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
ax.legend(loc=2,fontsize=5)
fig.set_size_inches(2.3,2)
#ax.set_ylim([1e-4,0.25])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([0.01,100])
ax.set_ylim([5e-4,10])
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
fig.savefig('Dpatch_newreac_controls_k2b_gef.pdf')
plt.show()


plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots() 
ax.errorbar(np.asarray(data6['params'])[[2,3]],Dp6,stdDp6,color='black',marker="o",markersize=2,capsize=2,label='GEF=200',linestyle='dotted')
ax.errorbar(data7['params'][0],Dp7,stdDp7,color='black',marker="o",markersize=2,capsize=2,label='GEF=400',linestyle='--')
ax.set_title(r'$k_{4a}=0$')
ax.set_xlabel(r'$k_{2b}$ $(s^{-1})$')
ax.set_xticks([1,3,10,30])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
ax.legend(loc=2,fontsize=5)
fig.set_size_inches(2.3,2)
#ax.set_ylim([1e-4,0.25])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([0.01,100])
ax.set_ylim([5e-4,10])
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
fig.savefig('Dpatch_newreac_controls_k4a_0_k2b_gef.pdf')
plt.show()


'''
#GAP rate vs Cdc42T , GAP rate vs Cdc42T-GEF
plt.close('all')
fig,ax = plt.subplots()


stddevtot42t= (np.asarray(data['mtot42t2']) - np.asarray(data['mtot42t'])**2)**0.5
stderrtot42t= stddevtot42t/((data['Nts']-data['skp'])*data['Ns'])**0.5
#ax.errorbar(data['mtot42t'],Dp,stdDp,stddevtot42t,color="black",linestyle='-',marker="o",markersize=2,capsize=2)

stddevtot42t2= (np.asarray(data2['mtot42t2']) - np.asarray(data2['mtot42t'])**2)**0.5
stderrtot42t2= stddevtot42t2/((data2['Nts']-data2['skp'])*data2['Ns'])**0.5
#ax.errorbar(data2['mtot42t'],Dp2,stdDp2,stddevtot42t2,color="red",linestyle='-',marker="o",markersize=2,capsize=2)

stddevtot42t3= (np.asarray(data3['mtot42t2']) - np.asarray(data3['mtot42t'])**2)**0.5
stderrtot42t3= stddevtot42t3/((data3['Nts']-data3['skp'])*data3['Ns'])**0.5
#ax.errorbar(data3['mtot42t'],Dp3,stdDp3,stddevtot42t3,color="blue",linestyle='-',marker="o",markersize=2,capsize=2)

stddevtot42t4= (np.asarray(data4['mtot42t2']) - np.asarray(data4['mtot42t'])**2)**0.5
stderrtot42t4= stddevtot42t4/((data4['Nts']-data4['skp'])*data4['Ns'])**0.5
#ax.errorbar(data4['mtot42t'],Dp4,stdDp4,stddevtot42t4,color="magenta",linestyle='-',marker="o",markersize=2,capsize=2)

ax.errorbar(data4['params'][0:3],data4['mtot42t'][0:3],stddevtot42t4[0:3],color='black',marker="o",markersize=2,capsize=2,label='GEF=50',linestyle='dotted')
ax.errorbar(data3['params'][0:3],data3['mtot42t'][0:3],stddevtot42t3[0:3],color='black',marker="o",markersize=2,capsize=2,label='GEF=100',linestyle='--')
#ax.errorbar(data2['params'],data2['mtot42t'],stddevtot42t2,color='red',marker="o",markersize=2,capsize=2,label='GEF=200')
ax.errorbar(data['params'],data['mtot42t'],stddevtot42t,color='black',marker="o",markersize=2,capsize=2,label='GEF=300',linestyle='-')

ax.set_xlabel(r'GAP rate $k_{2b}$ $(s^{-1})$')
ax.set_xticks([1,3,10,30])
ax.set_xlim([0.9,100])

ax.set_ylabel(r'Total Cdc42T')
#ax.set_ylim([0,0.25])
#plt.axhline(y=1./60,color='black',linestyle='--')
#ax.set_xlim([0.1,10])
ax.set_xscale('log')

ax2=ax.twinx()
# make a plot with different y-axis using second axis object


stddevgef42t= (np.asarray(data['mgef42t2']) - np.asarray(data['mgef42t'])**2)**0.5
stderrgef42t= stddevgef42t/((data['Nts']-data['skp'])*data['Ns'])**0.5
#ax.errorbar(data['mgef42t'],Dp,stdDp,stddevgef42t,color="black",linestyle='-',marker="o",markersize=2,capsize=2)

stddevgef42t2= (np.asarray(data2['mgef42t2']) - np.asarray(data2['mgef42t'])**2)**0.5
stderrgef42t2= stddevgef42t2/((data2['Nts']-data2['skp'])*data2['Ns'])**0.5
#ax.errorbar(data2['mgef42t'],Dp2,stdDp2,stddevgef42t2,color="red",linestyle='-',marker="o",markersize=2,capsize=2)

stddevgef42t3= (np.asarray(data3['mgef42t2']) - np.asarray(data3['mgef42t'])**2)**0.5
stderrgef42t3= stddevgef42t3/((data3['Nts']-data3['skp'])*data3['Ns'])**0.5
#ax.errorbar(data3['mgef42t'],Dp3,stdDp3,stddevgef42t3,color="blue",linestyle='-',marker="o",markersize=2,capsize=2)

stddevgef42t4= (np.asarray(data4['mgef42t2']) - np.asarray(data4['mgef42t'])**2)**0.5
stderrgef42t4= stddevgef42t4/((data4['Nts']-data4['skp'])*data4['Ns'])**0.5
#ax.errorbar(data4['mgef42t'],Dp4,stdDp4,stddevgef42t4,color="magenta",linestyle='-',marker="o",markersize=2,capsize=2)

ax2.errorbar(data4['params'][0:3],data4['mgef42t'][0:3],stddevgef42t4[0:3],color='red',marker="o",markersize=2,capsize=2,label='GEF=50',linestyle='dotted')
ax2.errorbar(data3['params'][0:3],data3['mgef42t'][0:3],stddevgef42t3[0:3],color='red',marker="o",markersize=2,capsize=2,label='GEF=100',linestyle='--')
#ax.errorbar(data2['params'],data2['mgef42t'],stddevgef42t2,color='red',marker="o",markersize=2,capsize=2,label='GEF=200')
ax2.errorbar(data['params'],data['mgef42t'],stddevgef42t,color='red',marker="o",markersize=2,capsize=2,label='GEF=300',linestyle='-')


ax2.set_ylabel("Cdc42T-GEF")
ax2.spines['right'].set_color('red')
ax2.tick_params(axis='y', colors='red')
#ax2.set_yticks([0,500,1000,1500,2000,2500])
ax2.yaxis.label.set_color('red')
#ax2.set_ylim([0,2700])
fig.set_size_inches(2.5,1.77)

plt.tight_layout()
fig.savefig('tot42t_gef42_gef_k2b.pdf')

'''

