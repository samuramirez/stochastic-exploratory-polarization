
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

#FRAP ANALYSIS


folder = './frap/'
N1=80
Nsims=30
T=90
samp=0.02

subfolders = ['dwellt42_newreac_gef','dwellt42_basegory_k2b_0p0625x_k5a_10x_gef']
labels=['High mobility','Low mobility']
parameters=[15,50,100,300,500 ] #gef
#fitintervalsizes=[3.6,5,6,6,6] #s for Cdc42m
fitintervalsizes=[[4.5,4.5,6,9,12],[30,30,30,30,30]] #for Cdc42T

#subfolders = ['frap42_basegory_k2b_0p0625x_k5a_10x_gef']
#parameters=[15,50,100,300,500 ] #gef
#fitintervalsizes=[1.5,1.5,1.5,1.5,1.5] #Cdc42m (s)
#fitintervalsizes=[30,30,30,30,30] #Cdc42T (s)



#eqconcinterval= np.arange(int(61.5/samp),int(62/samp)) #Cdc42m k2b/16 k5ax10
eqconcinterval= np.arange(int(80/samp),int(90/samp)) #Cdc42T k2b/16 k5ax10

stimulationt= int(60/samp)

meanmassinterval=np.arange(int(60/samp),int(90/samp))
 


col='black'

filenames = ['dynfrapN80gef'+str(parameters[i])+'sim' for i in range(len(parameters))]

#frap times
tausvsp=np.zeros([len(subfolders),len(parameters),3]) #mean, std, n
#mass
tot42tvsp=np.zeros([len(subfolders),len(parameters),3]) #mean, std, n
gef42tvsp=np.zeros([len(subfolders),len(parameters),3]) #mean, std, n
gefmvsp=np.zeros([len(subfolders),len(parameters),3]) #mean, std, n
cdc42tvsp=np.zeros([len(subfolders),len(parameters),3]) #mean, std, n
tot42vsp=np.zeros([len(subfolders),len(parameters),3]) #mean, std, n

#function to fit fluorescence decay
def frac_expl(x, tau):
    return np.exp(-x/tau)

#def frac_expl(x, tau1,tau2,f):
#    return 1 - f*np.exp(-x/tau1) - (1-f)*np.exp(-x/tau2)


#0 Cdc42T
#1 GEF42
#2 GEFm
#3 Cdc42m
#4 Cdc42c
#5 GEFc
#6 GEF42b
#7 Cdc42Tb
#8 Cdc42Db
#9 Cdc42Dcb

for p2,subfolder in enumerate(subfolders):
    
    fitinterval = [np.arange(int(60/samp),int((60+i)/samp)) for i in fitintervalsizes[p2]  ]
    times=[np.arange(len(fitinterval[i]))*samp for i in range(len(fitinterval)) ]
    TS_tot=[]
    for p,filename in enumerate(filenames):
        taus=[]
        n=0
        tot42tvst=np.zeros(int(T/samp))
        
        for sim in range(Nsims):
            print(sim)
            name =folder + subfolder +'/' + filename + str(sim+1) + '.dat'
          
            if os.path.isfile(name):
                #print(name)
                TS=np.genfromtxt(name)
                TS_tot.append(TS)
                #fluorescent Cdc42 at the membrane
                TS_gefm= np.sum(TS[:,(6,7)],1) #total Cdc42T 
                #TS_gefm= TS[:,3] #Cdc42Dm
                #subtract initial value at begining of fitinterval 
                #TS_gefm = TS_gefm - TS_gefm[fitinterval[p][0]]
                initconc=TS_gefm[stimulationt]
                
                #Get FRAP time
                    
                popt, pcov = curve_fit(frac_expl, times[p], TS_gefm[fitinterval[p]]/initconc)
                taus.append(popt[0])
                n+=1
                
                #mass
                tot42tvst = TS[:,0] + TS[:,1] + TS[:,6]
                tot42tvsp[p2,p,0]+= np.sum(tot42tvst[meanmassinterval])
                tot42tvsp[p2,p,1]+= np.sum(tot42tvst[meanmassinterval]**2)
                tot42tvsp[p2,p,2]+=len(meanmassinterval)
                
                #total membrane bound gef
                gefmvst = TS[:,1] + TS[:,2] + TS[:,6] + TS[:,7]
                #total gef42t
                #gefmvst = TS[:,1] + TS[:,6] 

                gefmvsp[p2,p,0]+= np.sum(gefmvst[meanmassinterval])
                gefmvsp[p2,p,1]+= np.sum(gefmvst[meanmassinterval]**2)
                gefmvsp[p2,p,2]+=len(meanmassinterval)
                              
            else:
                print "file " + name + " not found"
                
        tausvsp[p2,p,0]=np.mean(taus)
        tausvsp[p2,p,1]=np.std(taus)
        tausvsp[p2,p,2]=n
              



     
plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
#plt.errorbar(parameters,tausvsp[0,:,0], tausvsp[0,:,1],color='black',capsize=2)
fig,ax = plt.subplots()
ax.errorbar(parameters,tausvsp[1,:,0], 2*tausvsp[1,:,1]/tausvsp[1,:,2]**0.5,color='black',marker='s',markersize=2.5,capsize=2,label=labels[1])
ax.errorbar(parameters,tausvsp[0,:,0], 2*tausvsp[0,:,1]/tausvsp[0,:,2]**0.5,color='red',marker='.',markersize=2.5,capsize=2,label=labels[0])
ax.set_ylabel(r'Cdc42T dwell time (s)')
ax.set_xlabel(r'GEF')
ax.set_xlim([10,1000])
plt.xticks([10,100,200,300,400,500,1000])  
ax.set_ylim([1,300])
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend(loc=2,fontsize=6)    
fig.set_size_inches(2.3,1.9) 
plt.tight_layout()
fig.savefig('dwellt42_high_low_mobility_gef.pdf')

'''
plt.close('all')
matplotlib.rcParams.update({'font.size': 6})
#plt.errorbar(parameters,tausvsp[0,:,0], tausvsp[0,:,1],color='black',capsize=2)
fig2,ax2 = plt.subplots()
#ax.errorbar(parameters,tausvsp[0,:,0], tausvsp[0,:,1],color='black',capsize=2)
#ax2.set_xlabel(r'$k_{4a}$ ($\mu m^2/s$)')
ax2.set_xlim([5e-5,1])
ax2.set_xscale('log')
ax2.set_yscale('linear')
ax2.errorbar(parameters,gefmvsp[0,:,0]/gefmvsp[0,:,2], 2*(gefmvsp[0,:,1]/gefmvsp[0,:,2] - (gefmvsp[0,:,0]/gefmvsp[0,:,2])**2)**0.5/gefmvsp[0,:,2]**0.5,color='black',capsize=2)
ax2.set_ylabel(r'Total $GEF_{memb}$')
ax2.set_ylim([10,100])
ax2.set_yscale('linear')
fig2.set_size_inches(1.5,1)
plt.tight_layout()
#fig2.savefig('totgefm_newreac_gef100_k4a.pdf')
'''

'''
plt.close('all')
#plt.errorbar(parameters,tausvsp[0,:,0], tausvsp[0,:,1],color='black',capsize=2)
fig,ax = plt.subplots()
ax.errorbar(parameters,tausvsp[0,:,0], tausvsp[0,:,1]/tausvsp[0,:,2]**0.5,color='black')
#ax.errorbar(parameters,tausvsp[0,:,0], tausvsp[0,:,1],color='black',capsize=2)
ax.set_ylabel(r'FRAP $\tau$(s)')
ax.set_xlabel(r'$k_{4a}$ ($\mu m^2/s$)')
ax.set_xlim([5e-5,20])
ax.set_ylim([0.1,2.2])
ax.set_xscale('log')
ax.set_yscale('linear')
ax2=ax.twinx()
ax2.errorbar(parameters,gefmvsp[0,:,0]/gefmvsp[0,:,2], (gefmvsp[0,:,1]/gefmvsp[0,:,2] - (gefmvsp[0,:,0]/gefmvsp[0,:,2])**2)**0.5/gefmvsp[0,:,2]**0.5,color='red',capsize=2)
ax2.set_ylabel("Total membrane-bound GEF")
ax2.spines['right'].set_color('red')
#ax2.set_yticks([0,500,1000,1500,2000,2500])
ax2.yaxis.label.set_color('red')
ax2.set_ylim([10,100])
ax2.tick_params(axis='y', colors='red')
ax2.set_yscale('linear')
fig.set_size_inches(3.3,2.5)
plt.tight_layout()
fig.savefig('frap_totgefm_newreac_gef100_k4a.pdf')
'''

#fig=plt.figure(3)
#plt.errorbar(tot42tvsp[0,:,0]/tot42tvsp[0,:,2],tausvsp[0,:,0],yerr = tausvsp[0,:,1],xerr=(tot42tvsp[0,:,1]/tot42tvsp[0,:,2] - (tot42tvsp[0,:,0]/tot42tvsp[0,:,2])**2)**0.5,color='black',capsize=2)
#plt.errorbar(tot42tvsp[1,:,0]/tot42tvsp[1,:,2],tausvsp[1,:,0],yerr=tausvsp[1,:,1], xerr= (tot42tvsp[1,:,1]/tot42tvsp[1,:,2] - (tot42tvsp[1,:,0]/tot42tvsp[1,:,2])**2)**0.5,color='red',capsize=2)

'''
plt.close('all')
matplotlib.rcParams.update({'font.size': 7})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1

fig=plt.figure(1)
plt.errorbar(gefmvsp[0,:,0]/gefmvsp[0,:,2],tausvsp[0,:,0],yerr = tausvsp[0,:,1], xerr= (gefmvsp[0,:,1]/gefmvsp[0,:,2] - (gefmvsp[0,:,0]/gefmvsp[0,:,2])**2)**0.5,color='black',capsize=2, label=r'$k_3 = 10^{-2}\mu m^2/s$')
plt.errorbar(gefmvsp[1,:,0]/gefmvsp[1,:,2],tausvsp[1,:,0],yerr = tausvsp[1,:,1], xerr= (gefmvsp[1,:,1]/gefmvsp[1,:,2] - (gefmvsp[1,:,0]/gefmvsp[1,:,2])**2)**0.5,color='red',capsize=2,label=r'$k_3 = 10^{-1}\mu m^2/s $')
plt.xlabel('Membrane-bound GEF')
plt.ylabel(r'FRAP $\tau$ (s)')
plt.legend(loc=4,fontsize=5)
plt.xticks([0,100,200,300])
fig.set_size_inches(2.3,1.77)
#plt.ylim([0,20])
#plt.xlim([0,30])
plt.tight_layout()
#fig.savefig('frap_tau_memgef_k2a_0_k4a_0_k3_gef.pdf')
'''







TS_tot=np.asarray(TS_tot)



#plt.close('all')
fig=plt.figure(3)
#mass conservation
#0 Cdc42T
#1 GEF42
#2 GEFm
#3 Cdc42m
#4 Cdc42c
#5 GEFc
#6 GEF42b
#7 Cdc42Tb
#8 Cdc42Db
#9 Cdc42Dcb
#p=1
#plt.plot(np.sum(TS_tot[0,:,(0,1,3,4,6,7,8,9)],0))
#plt.plot(np.sum(TS_tot[p*30+0,:,(0,1,6,7)],0))
#plt.plot(TS_tot[p*30+2,:,1])
#plt.plot(np.sum(TS_tot[0,:,(1,2,5,6)],0))
#plt.plot(np.sum(TS_tot[0,:,(6,7,8,9)],0)) #marked molecules decay as they go to the cytosol and get unmarked
#plt.plot(TS_tot[0,:,9]) #no marked cytoplasmic
#plt.plot(np.sum(TS_tot[0,:,(0,1)],0)) #42t falls to almost zero as it is marked
#plt.plot(TS_tot[0,:,3]) #42dm falls to almost zero as it is marked
#plt.plot(TS_tot[0,:,8]) #marked Cdc42dm rapidly goes goes to the cytosol
#plt.plot(30*TS_tot[0,:,6]) #marked gef42 molecules decay as they go to the cytosol and get unmarked
#plt.plot(TS_tot[0,:,7]) #marked 42T molecules decay as they go to the cytosol and get unmarked


taus=[]
p=4
for i in range(30*p,30*p+5):
    #GEF membrane 
    TS_gefm= np.sum(TS_tot[i,:,(6,7,8)],0)
    #TS_gefm= TS_tot[i,:,3] #Cdc42m
    #TS_gefm = TS_gefm - TS_gefm[fitinterval[p][0]]
    initconc=TS_gefm[stimulationt]
    #eqconc=np.mean(TS_gefm[eqconcinterval])
    #popt, pcov = curve_fit(frac_expl, times[p], TS_gefm[fitinterval[p]]/eqconc,bounds=[[0,5,0],[1,20,0.8]])
    popt, pcov = curve_fit(frac_expl, times[p], TS_gefm[fitinterval[p]]/initconc)
    plt.plot(times[p], TS_gefm[fitinterval[p]]/initconc,color='black',alpha=0.5)
    plt.plot(times[p],frac_expl(times[p], *popt), '--',color='black')
    #taus.append(popt[0])
    #plt.plot(np.sum(TS_tot[i,:,(0,1)],0))

  
#GEF cytosol
#plt.plot(TS_tot[0,:,5])
#plt.ylim(0,1)
#plt.xlim(0,1)
#GEF bleached membrane 
#plt.plot(np.sum(TS_tot[0,:,(6,7)],0))
#GEF bleached cytosol
#plt.plot(TS_tot[0,:,8])



fig=plt.figure(4)
p=4
for i in range(30*p,30*p+5):
    #Cdc42 membrane 
    TS= np.sum(TS_tot[i,:,(6,7,8)],0)
    #TS= TS_tot[i,:,3]
    #Cdc42T GEF42
    #TS= np.sum(TS_tot[i,:,(0,1)],0)
    plt.plot(TS[int(0/samp):int(90/samp)],color='red',alpha=0.5)
    #plt.plot(TS)
    #plt.ylim(0,30)
