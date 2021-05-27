import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io
import os.path
from matplotlib import cm

plt.ion()
matplotlib.rcParams.update({'font.size': 12})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1.5

Nsims=30;
GEFconc=["200","300","400","500","600","700"];

#GEFconc=["500","600","700"];
#timeanalysis=[1:Nf];%verify this with H function file
#for i =1:length(GEFconc)

#Hallsims(l,:,:)=H_vals; %nsims x 50Hvalues * 30timepoints 

#verify this with matlab computation of H
max_r  = 5
dr=0.1
#in matlab r_vals = 0:dr:max_r;
r_vals = np.arange(0,max_r+dr,dr)
times=np.arange(10,310,10)
L=8
protr=0.02

#FOR FIGURE 4A
#meth=['kh','ksrhet'];
#Nmeth=["80"];
#figurename='H_PB_SG_h5rho_timeseries.pdf'
#cols= [ cm.jet(0.3) ,'red']

#FOR FIGURE 4C
meth=['ksrhet'];
Nmeth=["160"];
figurename='H_PB_SG_h2p5rho_timeseries.pdf'
cols= ['red']


methodname=[r'$k_h$',r'$k_c$']


#colors = [cm.jet(i) for i in np.linspace(1, 0, len(GEFconc))]
#colNs= [cm.jet(i) for i in np.linspace(1, 0, len(Nmeth))] 
#markers=['x','^']
#linestyles=['dotted','dashed','solid']
#linestyles=['dashed','solid']
gefconc='700'   

fi=1
plt.close('all')
# Clustering as a function of time: H(r==1.1) vs t
fig=plt.figure(fi)
fi+=1
#PARTICLE BASED
#Maximum of H function in particle-based simulations is consistently r=1.1
file =  filename='./HvaluesPB/Hvaluesgef' + gefconc + 'PB.mat';
if os.path.isfile(file):
    matdata = scipy.io.loadmat(file)
    
    Hallsims=matdata['Hallsims']
    #take average over multiple simulations and plot time series of H(r==1.)
    Hmeanvst=np.nanmean(Hallsims[:,(r_vals==1.1),:],0)[0]
    Hstdvst=np.nanstd(Hallsims[:,(r_vals==1.1),:],0)[0]
    plt.plot(times,Hmeanvst,color='black',label='particle-based',linewidth=2,linestyle='dashed')
    #plt.errorbar(times,Hmeanvst,Hstdvst,color='black',capsize=2)
    plt.fill_between(times,Hmeanvst-Hstdvst,Hmeanvst+Hstdvst,color='black',alpha=0.3)
    plt.xlabel(r'Time (s)')
    plt.ylabel(r'H(r=1.1$\mu m$)')
    #plt.ylim(2.5,3.5)
    #fig.set_size_inches(5, 4)
    #fig.set_size_inches(4, 3)
    #plt.legend(fontsize=8)
    plt.tight_layout()
else:
    print("file " + file + " not found")

#SPATIAL GILLESPIE    
for idxmet,met in enumerate(meth):
    for idxNmet,Nmet in enumerate(Nmeth):       
        file =  filename='./HvaluesSG/Hvalues'+'gef'+ gefconc +meth[idxmet]+'N'+Nmeth[idxNmet]+'.mat';       
        if os.path.isfile(file):
            matdata = scipy.io.loadmat(file)        
            Hallsims=matdata['Hallsims']      
            #look at average of individual simulations    
            #take average over multiple simulations and plot time series of H(r==1.)
            Hmeanvst=np.nanmean(Hallsims[:,(r_vals==1.1),:],0)[0]
            Hstdvst=np.nanstd(Hallsims[:,(r_vals==1.1),:],0)[0]
            plt.fill_between(times,Hmeanvst - Hstdvst,Hmeanvst + Hstdvst, color=cols[idxmet],alpha=0.3)
            plt.plot(times,Hmeanvst,color=cols[idxmet],label=methodname[idxmet]+ '  h = ' + '{:.1f}'.format(L/float(Nmeth[idxNmet])/protr)+r'$\rho$',linewidth=2)
            #plt.errorbar(times,Hmeanvst,Hstdvst,color=cols[idxmet],linestyle=linestyles[idxNmet],capsize=3)
            
        else:
            print("file " + file + " not found")
plt.legend(loc=4,fontsize=12)
plt.xlabel(r'Time (s)')
plt.ylabel(r'Clustering H(1.1$\mu m$)')
plt.ylim(0,3.8)
#plt.xlim(0,1.3)
fig.set_size_inches(4.2, 3.5)
plt.tight_layout()
#fig.savefig(figurename)
      
