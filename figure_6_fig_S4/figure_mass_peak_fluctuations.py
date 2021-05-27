import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io
import os.path
from matplotlib import cm
from statsmodels.graphics.tsaplots import plot_acf
from scipy.optimize import curve_fit
import pickle


folder='./'


subfolder = 'basegory_k2b_0p0625x_k5a_10x_gef_long'
with open(folder+subfolder+'_dataD.pkl', 'rb') as handle:
    datak2bk5a = pickle.load(handle)

folder = 'D:dynamic_polarity_data/'
subfolder = 'basegory_newreac_gef'
with open(folder+subfolder+'_dataD.pkl', 'rb') as handle:
    datak8 = pickle.load(handle)


folder='D:/dynamic_polarity_data/highres/'

recalculate=0 #1=calculate mass data. 0=load mass data, don't recalculate




subfolders = ['basegory_k2b_0p0625x_k5a_10x_gef/','basegory_newreac_gef/']


sims=[0,0]
variable='gef'
params=[15,25,100,300,500];
#p=0
colors = ['black','red','green']
labels=['Low mobility','High mobility','High mobility 2']
markers=['s','.']

Nspecies=2
T=900
sampling_sim = 1 #s sampling time in simulations
skipframes=0
nframes=T/sampling_sim - skipframes
times=(np.arange(nframes)+1)*sampling_sim/60.
L=8.
N1=80
Nsims=5
h=L/N1
taverpeak=np.arange(900-1*60,900)

if recalculate == 1:
    normcvpeak=np.zeros([len(subfolders),len(params),Nsims])
    totmass=np.zeros([len(subfolders),len(params),Nsims,nframes])
    totmassgef42=np.zeros([len(subfolders),len(params),Nsims,nframes])
    
    #taverpeak=np.arange(900-1*60,900)
    
    for sf,subfolder in enumerate(subfolders):
        for ip,p in enumerate(params):
            print(subfolder, p)
            for sim in range(Nsims):
                TS_cdc42T_gef42=0 #empty cdc42t distribution to avoid reusing when a .dat file is not present
                filename = 'dynN80'+variable+str(params[ip])+'sim'        
                name =folder + subfolder +'/' + filename + str(sim+1) + '.dat' 
                if os.path.isfile(name):
                    #print(name)
                    TS_cdc42T_gef42 = np.genfromtxt(name,skip_header=skipframes*Nspecies)
                    Nrows=np.shape(TS_cdc42T_gef42)[0]    
                    Ntimes=Nrows/Nspecies
                    if Ntimes != nframes:
                        print "Incorrect number of frames"
                        print(str(Ntimes), str(nframes), str(sim))       
                    N=int(np.sqrt(np.shape(TS_cdc42T_gef42)[1]))
                    if N1 != N:
                        print "Incorrect input N"   
                    #Get totalCdc42T = Cdc42T+GEF42 
                    TS_tot42T = TS_cdc42T_gef42[0::2,:] + TS_cdc42T_gef42[1::2,:]
                    TS_gef42T = TS_cdc42T_gef42[1::2,:]   
                    
                    #sum over grid element
                    #Tot 42T            
                    tot42tvst = np.sum(TS_tot42T[:,:],1)
                    #GEF42
                    gef42tvst = np.sum(TS_gef42T[:,:],1)
                    
                    totmass[sf,ip,sim,:]=tot42tvst
                    totmassgef42[sf,ip,sim,:]=gef42tvst
                    
                    #indices time, space                  
                    #average over time range
                    meantot42T_dist=np.mean(TS_tot42T[taverpeak,:],0)
                    stdtot42T_dist=np.std(TS_tot42T[taverpeak,:],0)
                    coefvar=stdtot42T_dist/meantot42T_dist
                    coefvar[np.isnan(coefvar)]=0
                    #mean cv of grid elements weighted by mean mass
                    normcvpeak[sf,ip,sim]=np.sum(coefvar*meantot42T_dist)/np.sum(meantot42T_dist)
                    #coefvars.append(coefvar)
                    #meandists.append(meantot42T_dist)
                    #stddists.append(stdtot42T_dist)
                    #put in matrix form
                    #meantot42T_mat = np.reshape(meantot42T_dist,(N,N))
                    #coefvar_mat = np.reshape(coefvar,(N,N))
                                
                else:
                    print "file " + name + " not found"
    
    
    #SAVE HIGH RESOLUTION DATA 
    data = {'normcvpeak':normcvpeak, 'totmass':totmass, 'totmassgef42':totmassgef42}
    with open(folder+'peakmassdata.pkl', 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    #LOAD HIGH RESOLUTION DATA
    folder='./highres/'
    with open(folder+'peakmassdata.pkl', 'rb') as handle:
        data=pickle.load(handle)
    normcvpeak=data['normcvpeak']
    totmass=data['totmass']
    totmassgef42=data['totmassgef42']


peakprofiles = np.zeros([len(subfolders),len(taverpeak),N1])     
for sf,subfolder in enumerate(subfolders):
    for ip,p in enumerate([25]):
        for sim in [1]:
            TS_cdc42T_gef42=0 #empty cdc42t distribution to avoid reusing when a .dat file is not present
            filename = 'dynN80'+variable+str(p)+'sim'        
            name =folder + subfolder +'/' + filename + str(sim+1) + '.dat' 
            if os.path.isfile(name):
                #print(name)
                TS_cdc42T_gef42 = np.genfromtxt(name,skip_header=skipframes*Nspecies)
                Nrows=np.shape(TS_cdc42T_gef42)[0]    
                Ntimes=Nrows/Nspecies
                if Ntimes != nframes:
                    print "Incorrect number of frames"
                    print(str(Ntimes), str(nframes), str(sim))       
                N=int(np.sqrt(np.shape(TS_cdc42T_gef42)[1]))
                if N1 != N:
                    print "Incorrect input N"   
                #Get totalCdc42T = Cdc42T+GEF42 
                TS_tot42T = TS_cdc42T_gef42[0::2,:] + TS_cdc42T_gef42[1::2,:]
                TS_gef42T = TS_cdc42T_gef42[1::2,:]   
                
                #get coordinates of max value in the first frame
                mat_first = np.reshape(TS_tot42T[taverpeak[0],:],(N,N))
                im,jm = np.unravel_index(mat_first.argmax(), mat_first.shape)
                            
                TS_centered=np.zeros((len(taverpeak),N*N))
                for ii,i in enumerate(taverpeak):
                    #set as matrix
                    mat=np.reshape(TS_tot42T[i,:],(N,N))
                    #center
                    centered=mat[np.ix_((np.arange(N) + im - int(N/2)) % N , (np.arange(N) + jm - int(N/2)) % N)]
                    TS_centered[ii,:]=centered.flat
                
                #meantot42T=np.mean(TS_centered[:,:],0)/h**2
                #meantot42T_mat.append(np.reshape(meantot42T,(N,N)))
                #get lateral profile of peaks
                for i in range(len(taverpeak)):
                    TS_centered_mat=np.reshape(TS_centered[i,:],(N,N))
                    peakprofiles[sf,i,:]=np.max(TS_centered_mat,0)
                         
            else:
                print "file " + name + " not found"
                
  

def estimated_autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result

def acf(series):
    n = len(series)
    data = np.asarray(series)
    mean = np.mean(data)
    c0 = np.sum((data - mean) ** 2) / float(n)

    def r(h):
        acf_lag = ((data[:n - h] - mean) * (data[h:] - mean)).sum() / float(n) / c0
        return round(acf_lag, 3)
    x = np.arange(n) # Avoiding lag 0 calculation
    acf_coeffs = map(r, x)
    return acf_coeffs


taver=np.arange(900-5*60,900)
intaver=[5,10,10,10,10]
tavers=[np.arange(900-i*60,900) for i in intaver]

meantotmass=np.zeros([len(subfolders),len(params),Nsims])
stdtotmass=np.zeros([len(subfolders),len(params),Nsims])
cvtotmass=np.zeros([len(subfolders),len(params),Nsims])

meantotmassgef42=np.zeros([len(subfolders),len(params),Nsims])
stdtotmassgef42=np.zeros([len(subfolders),len(params),Nsims])
cvtotmassgef42=np.zeros([len(subfolders),len(params),Nsims])

for sf,subfolder in enumerate(subfolders):
    for ip,p in enumerate(params):
        for sim in range(Nsims):
            meantotmass[sf,ip,sim]=np.mean(totmass[sf,ip,sim,tavers[ip]])
            stdtotmass[sf,ip,sim]=np.std(totmass[sf,ip,sim,tavers[ip]])
            cvtotmass[sf,ip,sim]=stdtotmass[sf,ip,sim]/meantotmass[sf,ip,sim]
            meantotmassgef42[sf,ip,sim]=np.mean(totmassgef42[sf,ip,sim,tavers[ip]])
            stdtotmassgef42[sf,ip,sim]=np.std(totmassgef42[sf,ip,sim,tavers[ip]])
            cvtotmassgef42[sf,ip,sim]=stdtotmassgef42[sf,ip,sim]/meantotmassgef42[sf,ip,sim]


plt.ion()
matplotlib.rcParams.update({'font.size': 7})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1


#AUTOCORRELATION CALCULATIONS
#function to fit
def expdecay(x, tau):
    return np.exp(-x/tau)
#k2b k5 points autocorrelation
npoints=np.zeros([len(subfolders),len(params)])
npoints[0,:]=[10, 10, 10, 10, 10  ]
#k8 points autocorrelation
npoints[1,:]=[ 4,5, 8, 8 ,8  ]
#k7 points autocorrelation
#npoints[2,:]=[10, 5, 5, 9, 7   ]
#plot_acf(totmass[1,4,1,tavers[2]])

taus=np.zeros([len(subfolders),len(params),Nsims])
for sf,subfolder in enumerate(subfolders):
    for ip,p in enumerate(params):
        for sim in range(Nsims):
            autocor=np.asarray(acf(totmass[sf,ip,sim,tavers[ip]]))
            times = np.arange(npoints[sf,ip])
            popt, pcov = curve_fit(expdecay, times, autocor[range(int(npoints[sf,ip]))])
            taus[sf,ip,sim]=popt[0]

sf=0
ip=3        
for sim in range(Nsims):
    autocor=np.asarray(acf(totmass[sf,ip,sim,taver]))
    times = np.arange(npoints[sf,ip])
    timesext= np.arange(npoints[sf,ip]+10)
    popt, pcov = curve_fit(expdecay, times, autocor[range(int(npoints[sf,ip]))])
    print(popt,pcov)
    plt.plot(timesext,autocor[ range(int(npoints[sf,ip]+10)) ] )
    plt.plot(timesext,expdecay(timesext, *popt), '--',alpha=0.5)
    plt.ylim([0,1])

sf=1
ip=4 
for sim in range(Nsims):
    plt.plot(totmass[sf,ip,sim,0:50])

fig2=plt.figure(2)
for sf in range(len(subfolders)-1):
    plt.errorbar(params,np.mean(meantotmass[sf,:,:],1),np.std(meantotmass[sf,:,:],1) ,linewidth=1,color=colors[sf],label=labels[sf],marker=markers[sf],markersize=3,capsize=2)
plt.ylim([0,5000])
plt.legend(loc=4,fontsize=6)    
plt.xlabel(r'GEF')
plt.ylabel(r'Mean total Cdc42T')
plt.xlim(10,1000)  
plt.xscale('log')
plt.xticks([10,100,200,300,400,500,1000])  
fig2.set_size_inches(2.3,1.9) 
fig2.tight_layout()
fig2.savefig('meantot42t_vs_GEF_stable_dynamic.pdf')    

fig10=plt.figure(10)
markers=['s','.']
for sf in range(len(subfolders)-1):
    plt.errorbar(params,np.mean(meantotmassgef42[sf,:,:],1),np.std(meantotmassgef42[sf,:,:],1) ,linewidth=1,color=colors[sf],label=labels[sf],marker=markers[sf],markersize=2.5,capsize=2)
#plt.ylim([0,5000])
plt.legend(loc=2,fontsize=6)    
plt.xlabel(r'GEF')
plt.ylabel(r'Mean Cdc42T-GEF')
plt.xlim(10,1000)  
plt.xscale('log')
plt.xticks([10,100,200,300,400,500,1000])  
fig10.set_size_inches(2.3,1.9) 
fig10.tight_layout()
fig10.savefig('meangef42t_vs_GEF_stable_dynamic.pdf')    



fig3=plt.figure(3)
plt.errorbar(params,np.mean(cvtotmass[0,:,1:5],1),np.std(cvtotmass[0,:,1:5],1) ,linewidth=1,color=colors[0],label=labels[0],marker=markers[0],markersize=2.5,capsize=2)
plt.errorbar(params,np.mean(cvtotmass[1,:,:],1),np.std(cvtotmass[1,:,:],1) ,linewidth=1,color=colors[1],label=labels[1],marker=markers[1],markersize=2.5,capsize=2)
plt.ylim([5e-3,1e-1])
#plt.yticks([0.005,1e-2,3e-2,1e-1])  
plt.yscale('log')
plt.legend(loc=1,fontsize=6)    
plt.xlabel(r'GEF')
plt.ylabel(r'CV total Cdc42T')
plt.xlim(10,1000)  
plt.xticks([10,100,200,300,400,500,1000])  
plt.xscale('log')
fig3.set_size_inches(2.3,1.9) 
fig3.tight_layout()
fig3.savefig('cvtot42t_vs_GEF_stable_dynamic.pdf') 

fig11=plt.figure(11)
plt.errorbar(params,np.mean(cvtotmassgef42[0,:,1:5],1),np.std(cvtotmassgef42[0,:,1:5],1) ,linewidth=1,color=colors[0],label=labels[0],marker=markers[0],markersize=2.5,capsize=2)
plt.errorbar(params,np.mean(cvtotmassgef42[1,:,:],1),np.std(cvtotmassgef42[1,:,:],1) ,linewidth=1,color=colors[1],label=labels[1],marker=markers[1],markersize=2.5,capsize=2)
#plt.ylim([5e-3,1e-1])
#plt.yticks([0.005,1e-2,3e-2,1e-1])  
plt.yscale('log')
plt.legend(loc=1,fontsize=6)    
plt.xlabel(r'GEF')
plt.ylabel(r'CV Cdc42T-GEF')
plt.xlim(10,1000)  
plt.xticks([10,100,200,300,400,500,1000])  
plt.xscale('log')
fig11.set_size_inches(2.3,1.9) 
fig11.tight_layout()
fig11.savefig('cvgef24_vs_GEF_stable_dynamic.pdf') 


plt.close('all')
fig4=plt.figure(4)
plt.errorbar(params,np.mean(taus[0,:,1:5],1),np.std(taus[0,:,1:5],1) ,linewidth=1,color=colors[0],label=labels[0],marker="o",markersize=2,capsize=2)
plt.errorbar(params,np.mean(taus[1,:,:],1),np.std(taus[1,:,:],1) ,linewidth=1,color=colors[1],label=labels[1],marker="o",markersize=2,capsize=2)
plt.ylim([1,100])
#plt.yticks([0.005,1e-2,3e-2,1e-1])  
plt.yscale('log')
plt.legend(loc=1,fontsize=6)    
plt.xlabel(r'GEF')
plt.ylabel(r'Cdc42T autocorrelation time (s)')
plt.xlim(10,1000)  
plt.xticks([10,100,200,300,400,500,1000])  
plt.xscale('log')
fig4.set_size_inches(2.3,1.9) 
fig4.tight_layout()
fig4.savefig('autocortau_GEF_stable_dynamic.pdf') 

fig5=plt.figure(5)
for sf in range(len(subfolders)-1):
    plt.errorbar(params,np.mean(normcvpeak[sf,:,:],1),np.std(normcvpeak[sf,:,:],1) ,linewidth=1,color=colors[sf],label=labels[sf],marker=markers[sf],markersize=2.5,capsize=2)
#plt.ylim([0,5000])
plt.legend(loc=1,fontsize=6)    
plt.xlabel(r'GEF')
plt.ylabel(r'$CV_{patch}$')
plt.xlim(10,1000)  
plt.xscale('log')
plt.yscale('linear')
plt.xticks([10,100,200,300,400,500,1000])  
fig5.set_size_inches(2.3,1.9) 
fig5.tight_layout()
fig5.savefig('CVpatch_vs_GEF_stable_dynamic.pdf')   



'''
fig5=plt.figure(5)
for sf in range(len(subfolders)-1):
    plt.errorbar(params,np.mean(stdtotmass[sf,:,1:5],1),np.std(stdtotmass[sf,:,1:5],1) ,linewidth=1,color=colors[sf],label=labels[sf],marker="o",markersize=2,capsize=2)
#plt.ylim([0,5000])
plt.legend(loc=4,fontsize=6)    
plt.xlabel(r'GEF')
plt.ylabel(r'Std dev total Cdc42T')
plt.xlim(10,1000)  
plt.xscale('log')
plt.xticks([10,100,200,300,400,500,1000])  
fig5.set_size_inches(2.3,1.9) 
fig5.tight_layout()
'''

for sf in range(len(subfolders)):
    plt.plot(params,taus[sf,:],linewidth=1)
plt.ylim([0,20])


#plt.figure()
#for sf in range(len(subfolders)):
#    plt.plot(params,cvtotmass[sf,:],linewidth=1)
#    plt.scatter(params,cvtotmass[sf,:],linewidth=1)
#plt.figure()
#for sf in range(len(subfolders)):
#    plt.plot(params,normcvpeak[sf,:],linewidth=1)




fig1 = plt.figure(1)
ax=fig1.add_subplot(1,1,1)
ax.errorbar(datak2bk5a['params'],datak2bk5a['Dp'],datak2bk5a['stdDp'],color='black',marker=markers[0],markersize=2.5,capsize=2,label=r'Low mobility',linestyle='-')
ax.errorbar(datak8['params'],datak8['Dp'],datak8['stdDp'],color='red',marker=markers[1],markersize=2.5,capsize=2,label=r'High mobility',linestyle='-')
ax.set_xlabel(r'GEF')
ax.set_xticks([10,100,200,300,400,500])
ax.set_ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
ax.legend(loc=1,fontsize=6)
fig1.set_size_inches(2.3,1.9)
ax.set_ylim([5e-4,1])
#plt.axhline(y=1./60,color='black',linestyle='--')
ax.set_xlim([10,1000])
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
fig1.savefig('Dpatch_stable_dynamic_gef.pdf')



Ds=[]
normcvpeaks=[]
meantotmasses=[]
stdtotmasses=[]
cvtotmasses=[]

#extract data from this parameters to plot with CVs
#params=[15,25,100,300,500];


#datak2bk5['params']
#[15, 25, 50, 100, 300, 500]
indsk2bk5a=[0,1,3,4,5]
for i in indsk2bk5a:
    Ds.append(datak2bk5a['Dp'][i])

for i in [0,1,2,3,4]:
    normcvpeaks.append(normcvpeak[0,i])
    meantotmasses.append(meantotmass[0,i])
    stdtotmasses.append(stdtotmass[0,i])
    cvtotmasses.append(cvtotmass[0,i])

#datak8['params']
#[15, 25, 50, 100, 200, 300, 500]
indsk8=[0,1,3,5,6]
for i in indsk8:
    Ds.append(datak8['Dp'][i])

for i in [0,1,2,3,4]:
    normcvpeaks.append(normcvpeak[1,i])
    meantotmasses.append(meantotmass[1,i])
    stdtotmasses.append(stdtotmass[1,i])
    cvtotmasses.append(cvtotmass[1,i])

#datak7['params']
#[15, 25, 50, 100, 200, 500]
#indsk7= [0,1,3,5]
#indscvpeakk7=[0,1,2,4]    
#for i in indsk7:
#    Ds.append(datak7['Dp'][i])

#for i in indscvpeakk7:
#    normcvpeaks.append(normcvpeak[2,i])
#    meantotmasses.append(meantotmass[2,i])
#    stdtotmasses.append(stdtotmass[2,i])
#    cvtotmasses.append(cvtotmass[2,i])



    
#plt.scatter(normcvpeaks,Ds)
#plt.figure(7)
#plt.scatter(normcvpeaks[0:5],Ds[0:5])
#plt.scatter(normcvpeaks[5:10],Ds[5:10])
#plt.scatter(normcvpeaks[10:14],Ds[10:14])
#plt.xscale('log')    
#plt.yscale('log')    

plt.close('all')
fig6=plt.figure(6)
plt.errorbar(np.mean(normcvpeak[0,:,:],1),np.asarray(datak2bk5a['Dp'])[indsk2bk5a],xerr=np.std(normcvpeak[0,:,:],1),yerr=np.asarray(datak2bk5a['stdDp'])[indsk2bk5a] ,linewidth=1,color=colors[0],label=labels[0],marker=markers[0],markersize=2.5,capsize=2)
plt.errorbar(np.mean(normcvpeak[1,:,:],1),np.asarray(datak8['Dp'])[indsk8],xerr=np.std(normcvpeak[1,:,:],1),yerr=np.asarray(datak8['stdDp'])[indsk8] ,linewidth=1,color=colors[1],label=labels[1],marker=markers[1],markersize=2.5,capsize=2)
#plt.errorbar(np.mean(normcvpeak[2,indscvpeakk7,:],1),np.asarray(Dpk7)[indsk7],xerr=np.std(normcvpeak[2,indscvpeakk7,:],1),yerr=np.asarray(stdDpk7)[indsk7] ,linewidth=1,color=colors[2],label=labels[2],marker="o",markersize=2,capsize=2)
plt.xscale('linear')
plt.yscale('log')
#plt.ylim([0,5000])
plt.legend(loc=4,fontsize=6)    
plt.xlabel(r'$CV_{patch}$')
plt.ylabel(r'$D_{patch}$ ($\mu m^2/min$)')
#plt.xlim(0.1,5)  
#plt.xticks([10,100,200,300,400,500,1000])  
fig6.set_size_inches(2.3,1.9) 
fig6.tight_layout()
fig6.savefig('Dp_vs_CVpatch_stable_dynamic.pdf')   


'''
plt.figure(2)
plt.scatter(meantotmasses[0:5],Ds[0:5])
plt.scatter(meantotmasses[5:10],Ds[5:10])
plt.scatter(meantotmasses[10:14],Ds[10:14])
plt.yscale('log')
plt.xscale('linear')    

plt.figure(3)
plt.scatter(stdtotmasses[0:5],Ds[0:5])
plt.scatter(stdtotmasses[5:10],Ds[5:10])
plt.scatter(stdtotmasses[10:14],Ds[10:14])
plt.yscale('log')
plt.xscale('log')   
 
plt.figure(4)
plt.scatter(cvtotmasses[0:5],Ds[0:5])
plt.scatter(cvtotmasses[5:10],Ds[5:10])
plt.scatter(cvtotmasses[10:14],Ds[10:14])
plt.yscale('log')
plt.xscale('log')    
'''



'''
plt.close('all')

fig1,ax1 = plt.subplots()
ax1.plot(times[0::2],totmass[0][0::2],color=colors[0],linewidth=0.5,label=labels[0])
ax1.plot(times[0::2],totmass[1][0::2],color=colors[1],linewidth=0.5,label=labels[1])   
ax1.legend(loc=4,fontsize=7)    
ax1.set_xlabel(r'Time (min)')
ax1.set_ylabel(r'Total mass')
#ax1.set_ylim(0,700)
#ax.set_ylim(0,3.5)
ax1.set_xlim(0,15)    
fig1.set_size_inches(2.3,1.9) 
fig1.tight_layout()
fig1.savefig('tot42t_stable_dynamic.pdf')

'''
plt.close('all')    
timesprofiles=[0,20,40,59]
fig9=plt.figure(9)
axes=[]
for ifig,i in enumerate(timesprofiles):
    axes.append(fig9.add_subplot(1,4,ifig+1))
    axes[ifig].plot(np.arange(0,8,0.1),peakprofiles[0,i,:]/h**2,color=colors[0])
    axes[ifig].set_xticks([0,4,8])
    axes[ifig].set_xlim(0,8)
    axes[ifig].set_xlabel(r'position ($\mu m$)')
    axes[ifig].set_ylabel(r'[Cdc42T] ($\mu m^{-2}$)')
    axes[ifig].set_ylim(0,1300)
    if ifig > 0:
        axes[ifig].set_yticks([])
        y_axis=axes[ifig].axes.get_yaxis()
        y_label = y_axis.get_label()
        y_label.set_visible(False)
fig9.set_size_inches(4.2,1.3)
fig9.tight_layout()
fig9.savefig('profiles_stable.pdf')

fig10=plt.figure(10)
axes=[]
for ifig,i in enumerate(timesprofiles):
    axes.append(fig10.add_subplot(1,4,ifig+1))
    axes[ifig].plot(np.arange(0,8,0.1),peakprofiles[1,i,:]/h**2,color=colors[1])
    axes[ifig].set_xticks([0,4,8])
    axes[ifig].set_ylim(0,8)
    axes[ifig].set_xlabel(r'position ($\mu m$)')
    axes[ifig].set_ylabel(r'[Cdc42T] ($\mu m^{-2}$)')
    axes[ifig].set_ylim(0,3500)
    axes[ifig].set_yticks([0,1000,2000,3000])

    if ifig > 0:
        axes[ifig].set_yticks([])
        y_axis=axes[ifig].axes.get_yaxis()
        y_label = y_axis.get_label()
        y_label.set_visible(False)
fig10.set_size_inches(4.2,1.3)
fig10.tight_layout()
fig10.savefig('profiles_dynamic.pdf')


