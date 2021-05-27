
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



folder='./'


subfolder = 'basegory_k7_1_gef'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data = pickle.load(handle)
    


folder='D:/dynamic_polarity_data/basegory/'

subfolder = 'basegory_k7_50_gef_long'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data4 = pickle.load(handle)

subfolder = 'basegory_k7_10_gef_long'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data2 = pickle.load(handle)

subfolder = 'basegory_k7_20_gef_long'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data3 = pickle.load(handle)



    
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


#MSD k7=1
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp=[]
a=[]
stdDp=[]
stda=[]
nlinear=12 #?
n0=1   
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
    plt.xlim([0,20])
    plt.yscale('log')
    plt.xscale('log')
    #plt.legend(loc=2,fontsize=12)
    fig.set_size_inches(2.5,1.77)   
#plt.title('Max step considered='+str(maxstep))
plt.tight_layout()


#MSD k7=10
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp2=[]
a2=[]
stdDp2=[]
stda2=[]
for i in range(len(data2['params'])):
    if data2['params'][i] == 500  :
        n0=4
        nlinear=12
    else:
        n0=1
        nlinear=7
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

#MSD k7=20
plt.close('all')
fi=1  
fig=plt.figure(fi)
fi+=1
Dp3=[]
a3=[]
stdDp3=[]
stda3=[]
for i in range(len(data3['params'])):
    if data3['params'][i] == 500  :
        n0=4
        nlinear=15
    else:
        n0=1
        nlinear=8
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



#MSD k7 = 50
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
    
    if data4['params'][i] == 500 or data4['params'][i] == 200 :
        n0=6
        nlinear=10
    else:
        n0=1
        nlinear=7
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


#THIS IS FIGURE 5E.  Dpatch vs GEF for k4a=2 and k4a=0 LOG
plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots() 
#ax.plot(parameters,Dp,color='black',marker="o",markersize=3)
ax.errorbar(data['params'],Dp,stdDp,color='black',marker="o",markersize=2,capsize=2,label=r'$k_{6}=1\mu m^2/s$ ',linestyle='-')
#ax.errorbar(data2['params'],Dp2,stdDp2,color='black',marker="o",markersize=2,capsize=2,label=r'$k_{6}=10\mu m^2/s$ ',linestyle='-')
#ax.errorbar(data3['params'],Dp3,stdDp3,color='black',marker="o",markersize=2,capsize=2,label=r'$k_{6}=20\mu m^2/s$ ')
ax.errorbar(data4['params'],Dp4,stdDp4,color='black',marker="o",markersize=2,capsize=2,label=r'$k_{6}=50\mu m^2/s$ ',linestyle='-')

ax.set_xlabel(r'GEF')
ax.set_xticks([0,100,200,300,400,500])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')

#ax2=ax.twinx()
# make a plot with different y-axis using second axis object
#stddevgef42t= (np.asarray(data4['mgef42t2']) - np.asarray(data4['mgef42t'])**2)**0.5
#ax2.errorbar(data4['params'],data4['mgef42t'],stddevgef42t,color='orange',marker="o",markersize=2,capsize=2,linestyle='-')

#ax.legend(loc=1,fontsize=7)
ax.set_ylim([5e-4,1])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([0,510])
#ax.set_xscale('log')
ax.set_yscale('log')

#ax2.set_ylabel("Cdc42T-GEF")
#ax2.spines['right'].set_color('orange')
#ax2.tick_params(axis='y', colors='orange')
#ax2.set_yticks([0,500,1000,1500,2000,2500])
#ax2.yaxis.label.set_color('orange')
#ax2.set_ylim([0,500])

fig.set_size_inches(2.3,1.77)
#fig.set_size_inches(2.5,1.77) # after adding second axis

plt.tight_layout()
#fig.savefig('Dpatch_k7_1_50_gef.pdf')
plt.show()

#SAVE DATA
subfolder = 'basegory_k7_50_gef_long'
dataD = {'params':data4['params'], 'Dp':Dp4, 'stdDp':stdDp4, 'a':a4, 'stda':stda4}
with open(folder+subfolder+'_dataD.pkl', 'wb') as handle:
    pickle.dump(dataD, handle, protocol=pickle.HIGHEST_PROTOCOL)




'''
#Dpatch vs Cdc42T-GEF
plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots()
#stddevtot42t= (np.asarray(data['mtot42t2']) - np.asarray(data['mtot42t'])**2)**0.5
#stderrtot42t= stddevtot42t/((data['Nts']-data['skp'])*data['Ns'])**0.5
#ax.errorbar(data['mtot42t'],Dp,stdDp,color='black',marker="o",markersize=2,capsize=2,label=r'$k_{4a}=2$, $k_{7}=0.5\mu m^2/s$ ')
#stddevtot42t2= (np.asarray(data2['mtot42t2']) - np.asarray(data2['mtot42t'])**2)**0.5
#stderrtot42t2= stddevtot42t2/((data2['Nts']-data2['skp'])*data2['Ns'])**0.5
#ax.errorbar( data2['mtot42t'],Dp2,stdDp2,color="red",linestyle='-',marker="o",markersize=2,capsize=2)
#ax.set_xlabel(r'Total Cdc42T')
#ax.set_xlim([0,4000])


stddevgef42t= (np.asarray(data['mgef42t2']) - np.asarray(data['mgef42t'])**2)**0.5
stderrgef42t= stddevgef42t/((data['Nts']-data['skp'])*data['Ns'])**0.5
ax.errorbar(data['mgef42t'],Dp,stdDp,stddevgef42t,color="black",linestyle='-',marker="o",markersize=2,capsize=2)

stddevgef42t2= (np.asarray(data2['mgef42t2']) - np.asarray(data2['mgef42t'])**2)**0.5
stderrgef42t2= stddevgef42t2/((data2['Nts']-data2['skp'])*data2['Ns'])**0.5
ax.errorbar(data2['mgef42t'],Dp2,stdDp2,stddevgef42t2,color="red",linestyle='-',marker="o",markersize=2,capsize=2)

stddevgef42t3= (np.asarray(data3['mgef42t2']) - np.asarray(data3['mgef42t'])**2)**0.5
stderrgef42t3= stddevgef42t3/((data3['Nts']-data3['skp'])*data3['Ns'])**0.5
ax.errorbar(data3['mgef42t'],Dp3,stdDp3,stddevgef42t3,color="blue",linestyle='-',marker="o",markersize=2,capsize=2)

stddevgef42t4= (np.asarray(data4['mgef42t2']) - np.asarray(data4['mgef42t'])**2)**0.5
stderrgef42t4= stddevgef42t4/((data4['Nts']-data4['skp'])*data4['Ns'])**0.5
ax.errorbar(data4['mgef42t'],Dp4,stdDp4,stddevgef42t4,color="magenta",linestyle='-',marker="o",markersize=2,capsize=2)


ax.set_xlabel(r'Cdc42T-GEF')
ax.set_xlim([0,200])


ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
#ax.set_xticks([0,100,300,500])
#ax.legend(loc=1,fontsize=7)
fig.set_size_inches(2.3,1.77)
#ax.set_ylim([1e-4,0.25])
#plt.axhline(y=1./60,color='black',linestyle='--')
#ax.set_xscale('log')
ax.set_yscale('log')

plt.tight_layout()
#fig.savefig('Dpatch_newreac_mm_0_GEF42_log.pdf')
plt.show()


'''
