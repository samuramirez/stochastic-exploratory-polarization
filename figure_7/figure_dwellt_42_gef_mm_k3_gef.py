
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
Nsims=10
T=90
samp=0.02

#subfolders = ['dwellt42_newreac_k3_e-2_gef','dwellt42_newreac_k3_e-1_gef','dwelltgef_newreac_k3_e-2_gef','dwelltgef_newreac_k3_e-1_gef']
#xvalues=[1e-2,1e-1]
#xlabel=r'$k_{3}$ ($\mu m^2/s$)'
#rate='k3'


#subfolders = ['dwellt42_newreac_k2_e-2_gef','dwellt42_newreac_k2_e-1_gef','dwelltgef_newreac_k2_e-2_gef','dwelltgef_newreac_k2_e-1_gef']
#xvalues=[1e-2,1e-1]
#xlabel=r'$k_{2a}$ ($\mu m^2/s$)'
#rate='k2a'

subfolders = ['dwellt42_newreac_k4a_e-3_gef','dwellt42_newreac_k4a_e-2_gef','dwelltgef_newreac_k4a_e-3_gef','dwelltgef_newreac_k4a_e-2_gef']
xlabel=r'$k_{4a}$ ($\mu m^2/s$)'
xvalues=[1e-3,1e-2]
rate='k4a'

color2='royalblue'
labels=['GEF=100','GEF=200','GEF=300']
#subfolders = ['dwelltgef_newreac_k2_e-2_gef','dwelltgef_newreac_k2_e-1_gef']
#labels=['k3 = 10^-2','k3 = 10^-1']
parameters=[100,200,300,400 ] #gef


#fitintervalsizes=[3.6,5,6,6,6] #s for Cdc42m
#fitintervalsizes=[[4.5,4.5,6,9,12],[30,30,30,30,30]] #for Cdc42T
fitintervalsizes=[[4.5,4.5,6,6],[4.5,6,7.5,9],[1,1,1,1],[1,1,1,1]] #for Cdc42T

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

#function to fit fluoresence decay
def frac_expl(x, tau):
    return np.exp(-x/tau)

#def frac_expl(x, tau1,tau2,f):
#    return 1 - f*np.exp(-x/tau1) - (1-f)*np.exp(-x/tau2)

#Cdc42 dwellt sims
#0 Cdc42T
#1 GEF42
#2 GEFm
#3 Cdc42m
#4 Cdc42c
#5 GEFc
#6 GEF42b track this 
#7 Cdc42Tb track this
#8 Cdc42Db #is faster compared to gef42 and cdc42T
#9 Cdc42Dcb #is zero

#GEF dwellt sims
#0 Cdc42T
#1 GEF42
#2 GEFm
#3 Cdc42m
#4 Cdc42c
#5 GEFc
#6 GEF42b track this
#7 GEFmb  track this
#8 GEFcb #is zero



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
                TS_gefm= np.sum(TS[:,(6,7)],1) #total marked protein (works for gef and cdc42 dwellt sims)
                #TS_gefm= TS[:,3] #Cdc42Dm
                #subtract initial value at begining of fitinterval 
                #TS_gefm = TS_gefm - TS_gefm[fitinterval[p][0]]
                initconc=TS_gefm[stimulationt]
                if initconc != 0:
                #Get FRAP time                    
                    popt, pcov = curve_fit(frac_expl, times[p], TS_gefm[fitinterval[p]]/initconc)
                    taus.append(popt[0])
                    n+=1
                
                #total cdc42t
                tot42tvst = TS[:,0] + TS[:,1] + TS[:,6] + TS[:,7]
                tot42tvsp[p2,p,0]+= np.sum(tot42tvst[meanmassinterval])
                tot42tvsp[p2,p,1]+= np.sum(tot42tvst[meanmassinterval]**2)
                tot42tvsp[p2,p,2]+=len(meanmassinterval)
                
                #total membrane bound gef
                #gefmvst = TS[:,1] + TS[:,2] + TS[:,6] + TS[:,7]
                
                #total gef42t
                gefmvst = TS[:,1] + TS[:,6] 
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


#DWELL TIMES
fig1 = plt.figure(1)
ax=fig1.add_subplot(111)
#fig,ax = plt.subplots()
#CDC42
#tausvsp indices: folder (dwellt42 gef, k3), gef, statistic
#ax.errorbar(xvalues,tausvsp[(0,1),0,0], 2*tausvsp[(0,1),0,1]/tausvsp[(0,1),0,2]**0.5,color='black',capsize=2,label=labels[0],linestyle='dotted')
#ax.errorbar(xvalues,tausvsp[(0,1),1,0], 2*tausvsp[(0,1),1,1]/tausvsp[(0,1),1,2]**0.5,color='black',capsize=2,label=labels[1],linestyle='--')
#tot gef = 300
ax.errorbar(xvalues,tausvsp[(0,1),2,0], 2*tausvsp[(0,1),2,1]/tausvsp[(0,1),2,2]**0.5,color='black',capsize=2,label=labels[2],linestyle='-')
ax.set_xlabel(xlabel)
#ax.set_xlim([0,0.11])
plt.xticks(xvalues)
ax.set_ylabel(r'Cdc42T dwell time (s)')
ax.set_ylim([1,3])
#ax.set_xscale('log')
#ax.set_yscale('log')
#plt.legend(loc=2,fontsize=6)    

#GEF dwell time
ax2=ax.twinx()
#ax.errorbar(xvalues,tausvsp[(2,3),0,0], 2*tausvsp[(2,3),0,1]/tausvsp[(2,3),0,2]**0.5,color='red',capsize=2,label=labels[0],linestyle='dotted')
#ax.errorbar(xvalues,tausvsp[(2,3),1,0], 2*tausvsp[(2,3),1,1]/tausvsp[(2,3),1,2]**0.5,color='red',capsize=2,label=labels[1],linestyle='--')
ax2.errorbar(xvalues,tausvsp[(2,3),2,0], 2*tausvsp[(2,3),2,1]/tausvsp[(2,3),2,2]**0.5,color=color2,capsize=2,label=labels[2],linestyle='-')
ax2.set_ylabel(r'GEF dwell time (s)')
ax2.spines['right'].set_color(color2)
ax2.tick_params(axis='y', colors=color2)
#ax2.set_yticks([0,500,1000,1500,2000,2500])
ax2.yaxis.label.set_color(color2)
ax2.set_ylim([0,0.6])
fig1.set_size_inches(1.5,1.9) 
plt.tight_layout()
fig.savefig('dwellt_42t_gef42_gef_300_'+rate+'.pdf')




#MASS AT THE CLUSTER
fig,ax = plt.subplots()
#Cdc42T
#ax.errorbar(xvalues,tot42tvsp[(0,1),0,0]/tot42tvsp[(0,1),0,2], 2*(tot42tvsp[(0,1),0,1]/tot42tvsp[(0,1),0,2] - (tot42tvsp[(0,1),0,0]/tot42tvsp[(0,1),0,2])**2)**0.5/tot42tvsp[(0,1),0,2]**0.5,color='black',capsize=2,label=labels[0])
#ax.errorbar(xvalues,tot42tvsp[(0,1),1,0]/tot42tvsp[(0,1),1,2], 2*(tot42tvsp[(0,1),1,1]/tot42tvsp[(0,1),1,2] - (tot42tvsp[(0,1),1,0]/tot42tvsp[(0,1),1,2])**2)**0.5/tot42tvsp[(0,1),1,2]**0.5,color='black',capsize=2,label=labels[1])
ax.errorbar(xvalues,tot42tvsp[(0,1),2,0]/tot42tvsp[(0,1),2,2], 2*(tot42tvsp[(0,1),2,1]/tot42tvsp[(0,1),2,2] - (tot42tvsp[(0,),2,0]/tot42tvsp[(0,1),2,2])**2)**0.5/tot42tvsp[(0,1),2,2]**0.5,color='black',capsize=2,label=labels[2])
ax.set_ylabel(r'Cdc42T')
ax.set_xlabel(xlabel)
#ax.set_xlim([90,410])
plt.xticks(xvalues)  
ax.set_ylim([0,3000])
#ax.set_xscale('log')
#ax.set_yscale('log')
#GEF42
ax2=ax.twinx()
#ax.errorbar(xvalues,gefmvsp[(0,1),0,0]/gefmvsp[(0,1),0,2], 2*(gefmvsp[(0,1),0,1]/gefmvsp[(0,1),0,2] - (gefmvsp[(0,1),0,0]/gefmvsp[(0,1),0,2])**2)**0.5/gefmvsp[(0,1),0,2]**0.5,color='red',capsize=2,label=labels[0])
#ax.errorbar(xvalues,gefmvsp[(0,1),1,0]/gefmvsp[(0,1),1,2], 2*(gefmvsp[(0,1),1,1]/gefmvsp[(0,1),1,2] - (gefmvsp[(0,1),1,0]/gefmvsp[(0,1),1,2])**2)**0.5/gefmvsp[(0,1),1,2]**0.5,color='red',capsize=2,label=labels[1])
ax2.errorbar(xvalues,gefmvsp[(0,1),2,0]/gefmvsp[(0,1),2,2], 2*(gefmvsp[(0,1),2,1]/gefmvsp[(0,1),2,2] - (gefmvsp[(0,),2,0]/gefmvsp[(0,1),2,2])**2)**0.5/gefmvsp[(0,1),2,2]**0.5,color=color2,capsize=2,label=labels[2])
ax2.set_ylabel(r'Cdc42T-GEF')
ax2.spines['right'].set_color(color2)
ax2.tick_params(axis='y', colors=color2)
#ax2.set_yticks([0,500,1000,1500,2000,2500])
ax2.yaxis.label.set_color(color2)
ax2.set_ylim([0,200])
#plt.legend(loc=2,fontsize=6)    
fig.set_size_inches(1.7,1.9) 
plt.tight_layout()
fig.savefig('mass_42t_gef42_gef_300_'+rate+'.pdf')

