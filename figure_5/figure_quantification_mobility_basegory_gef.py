
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


subfolder = 'basegory_newreac0_gef'
with open(folder+subfolder+'.pkl', 'rb') as handle:
    data = pickle.load(handle,encoding='latin1')


    
#data = {'ssd':sum_squared_disps, 'nsq':n_squared_disps, 'params':parameters,
#        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
#       
#        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
#       }

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
    plt.xlim([0,30])
    plt.yscale('log')
    plt.xscale('log')
    #plt.legend(loc=2,fontsize=12)
    fig.set_size_inches(2.5,1.77)   
#plt.title('Max step considered='+str(maxstep))
plt.tight_layout()



#THIS IS FIGURE 5E.  Dpatch vs GEF for k4a=2 and k4a=0 LOG
plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
fig,ax = plt.subplots() 
#ax.plot(parameters,Dp,color='black',marker="o",markersize=3)
ax.errorbar(data['params'],Dp,stdDp,color='black',marker="o",markersize=2,capsize=2,label=r'$k_{6}=1\mu m^2/s$ ',linestyle='-')

#ax.errorbar(data3['params'],Dp3,stdDp3,color='black',marker="o",markersize=2,capsize=2,label=r'$k_{6}=100\mu m^2/s$ ')
#ax.errorbar(data4['params'],Dp4,stdDp4,color='red',marker="o",markersize=2,capsize=2,label=r'$0.0625x \: k_{2b}$',linestyle='--')
#ax.errorbar(data5['params'],Dp5,stdDp5,color='red',marker="o",markersize=2,capsize=2,label=r'$0.0625x \: k_{2b}$',)

ax.set_xlabel(r'GEF')
ax.set_xticks([0,200,400,600])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
#ax.legend(loc=1,fontsize=7)
fig.set_size_inches(2.3,1.77)
ax.set_ylim([0,0.005])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([0,720])
#ax.set_xscale('log')
ax.set_yscale('linear')
plt.tight_layout()
#fig.savefig('Dpatch_basegory_gef.pdf')
plt.show()


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
#stderrgef42t= stddevgef42t/((60-data['skp'])*data['Ns'])**0.5
ax.errorbar(data['mgef42t'],Dp,stdDp,stddevgef42t,color="black",linestyle='-',marker="o",markersize=2,capsize=2)

stddevgef42t2= (np.asarray(data2['mgef42t2']) - np.asarray(data2['mgef42t'])**2)**0.5
#stderrgef42t2= stddevgef42t2/((60-data2['skp'])*data2['Ns'])**0.5
ax.errorbar(data2['mgef42t'],Dp2,stdDp2,stddevgef42t2,color="red",linestyle='-',marker="o",markersize=2,capsize=2)


#stddevgef42t3= (np.asarray(data3['mgef42t2']) - np.asarray(data3['mgef42t'])**2)**0.5
#stderrgef42t3= stddevgef42t3/((data3['Nts']-data3['skp'])*data3['Ns'])**0.5
#ax.errorbar(data3['mgef42t'],Dp3,stdDp3,stddevgef42t3,color="blue",linestyle='-',marker="o",markersize=2,capsize=2)

#stddevgef42t4= (np.asarray(data4['mgef42t2']) - np.asarray(data4['mgef42t'])**2)**0.5
#stderrgef42t4= stddevgef42t4/((data4['Nts']-data4['skp'])*data4['Ns'])**0.5
#ax.errorbar(data4['mgef42t'],Dp4,stdDp4,stddevgef42t4,color="magenta",linestyle='-',marker="o",markersize=2,capsize=2)

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
