import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io
import os.path
from matplotlib import cm
import pickle

plt.ion()
matplotlib.rcParams.update({'font.size': 12})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1.5

Nsims=3;

variable='gef'
params=[50,100,300,500];

folders=['../basegory/basegory_k3_0p1_gef/','../basegory/basegory_k3_1_gef/','../basegory/basegory_k3_10_gef/']
labels=[r'$k_{3}=0.1 \mu m^2/s$',r'$k_{3}=1 \mu m^2/s$',r'$k_{3}=10 \mu m^2/s$']
figurename="H_k3_gef.pdf"

folders=['../basegory/basegory_k4a_1_gef/','../basegory/basegory_k4a_10_gef/']
labels=[r'$k_{4a}=1 \mu m^2/s$',r'$k_{4a}=10\mu m^2/s$']
figurename="H_k4a_gef.pdf"

folders=['../basegory/basegory_k5a_3x_gef/','../basegory/basegory_k5a_10x_gef/','../basegory/basegory_k5a_30x_gef/','../basegory/basegory_k5a_100x_gef/']
labels=[r'$3x \: k_{5a}$',r'$10x \: k_{5a}$',r'$30x \: k_{5a}$',r'$100x \: k_{5a}$']
figurenameH="H_k5a_gef.pdf"
figurenameMass="mass_k5a_gef.pdf"



folders=['../basegory/basegory_k2a_0p1_gef','../basegory/basegory_k2a_1_gef','../basegory/basegory_k2a_10_gef']
title=r'$k_{2a} (\mu m^2/s)$'
labels=[r'0.1','1','10']
figurenameH="H_k2a_gef.pdf"
figurenameMass="mass_k2a_gef.pdf"

folders=['../basegory/basegory_k1a_3x_gef','../basegory/basegory_k1a_10x_gef','../basegory/basegory_k1a_30x_gef','../basegory/basegory_k1a_100x_gef']
title=r'$k_{1a}$'
labels=['3x','10x','30x','100x']
figurenameH="H_k1a_gef.pdf"
figurenameMass="mass_k1a_gef.pdf"


folders=['../basegory/basegory_k2b_0p5x_gef','../basegory/basegory_k2b_0p25x_gef','../basegory/basegory_k2b_0p125x_gef','../basegory/basegory_k2b_0p0625x_gef']
title=r'$k_{2b}$'
labels=['0.5x','0.25x','0.125','0.0625x']
figurenameH="H_k2b_gef.pdf"
figurenameMass="mass_k2b_gef.pdf"

folders=['../basegory/basegory_k3_0p1_gef','../basegory/basegory_k3_1_gef','../basegory/basegory_k3_10_gef']
title=r'$k_{3} (\mu m^2/s)$'
labels=['0.1','1','10']
figurenameH="H_k3_gef.pdf"
figurenameMass="mass_k3_gef.pdf"

folders=['../basegory/basegory_k4a_1_gef','../basegory/basegory_k4a_10_gef']
title=r'$k_{4a} (\mu m^2/s)$'
labels=['1','10']
figurenameH="H_k4a_gef.pdf"
figurenameMass="mass_k4a_gef.pdf"


folders=['../basegory/basegory_k5a_3x_gef','../basegory/basegory_k5a_10x_gef','../basegory/basegory_k5a_30x_gef','../basegory/basegory_k5a_100x_gef']
title=r'$k_{5a}$'
labels=['3x','10x','30x','100x']
figurenameH="H_k5a_gef.pdf"
figurenameMass="mass_k5a_gef.pdf"

folders=['D:dynamic_polarity_data/basegory/basegory_k7_1_gef','D:dynamic_polarity_data/basegory/basegory_k7_10_gef','D:dynamic_polarity_data/basegory/basegory_k7_20_gef','D:dynamic_polarity_data/basegory/basegory_k7_50_gef']
title=r'$k_{6}$'
labels=['5x','50x','100x','200x']
figurenameH="H_k7_gef.pdf"
figurenameMass="mass_k7_gef.pdf"
params=[15,25,50,100,300,500];

folders=['D:dynamic_polarity_data/basegory/basegory_k1b_0p5x_gef','D:dynamic_polarity_data/basegory/basegory_k1b_0p25x_gef','D:dynamic_polarity_data/basegory/basegory_k1b_0p125x_gef','D:dynamic_polarity_data/basegory/basegory_k1b_0p0625x_gef']
title=r'$k_{1b}$'
labels=['0.5x','0.25x','0.125','0.0625x']
figurenameH="H_k1b_gef.pdf"
figurenameMass="mass_k1b_gef.pdf"
params=[50,100,300,500];

folders=['D:dynamic_polarity_data/basegory/basegory_k5b_0p5x_gef','D:dynamic_polarity_data/basegory/basegory_k5b_0p25x_gef','D:dynamic_polarity_data/basegory/basegory_k5b_0p125x_gef','D:dynamic_polarity_data/basegory/basegory_k5b_0p0625x_gef']
title=r'$k_{5b}$'
labels=['0.5x','0.25x','0.125','0.0625x']
figurenameH="H_k5b_gef.pdf"
figurenameMass="mass_k5b_gef.pdf"
params=[50,100,300,500];

folders=['D:dynamic_polarity_data/basegory/basegory_k4b_0p5x_gef','D:dynamic_polarity_data/basegory/basegory_k4b_0p25x_gef','D:dynamic_polarity_data/basegory/basegory_k4b_0p125x_gef','D:dynamic_polarity_data/basegory/basegory_k4b_0p0625x_gef']
title=r'$k_{4b}$'
labels=['0.5x','0.25x','0.125','0.0625x']
figurenameH="H_k4b_gef.pdf"
figurenameMass="mass_k4b_gef.pdf"
params=[50,100,300,500];

meanHs=np.zeros([len(folders),len(params)])
stdHs=np.zeros([len(folders),len(params)])
meanmolec=np.zeros([len(folders),len(params)])
stddevmolec=np.zeros([len(folders),len(params)])
stderrmolec=np.zeros([len(folders),len(params)])


#timeanalysis=[1:Nf];%verify this with H function file

#Hallsims(l,:,:)=H_vals; %nsims x 50Hvalues * 30timepoints 

#verify this with matlab computation of H
max_r  = 5
dr=0.1
#in matlab r_vals = 0:dr:max_r;
r_vals = np.arange(0,max_r+dr,dr)
times=np.arange(1,21)
L=8
#10 min sims sampling 30s
skip=10 #to give time for polarity patch to form

#Maximum of H function consistently r=1.1
for f,folder in enumerate(folders):
    for p,param in enumerate(params):
        
        file =  filename=folder+'/Hvalues'+variable +str(param)+'N80.mat';
        if os.path.isfile(file):
            matdata = scipy.io.loadmat(file)
            
            Hallsims=matdata['Hallsims']
            #take average over multiple simulations multiple time points of H(r==1.)
            meanHs[f,p]=np.nanmean(Hallsims[:,(r_vals==1.1),skip:])
            stdHs[f,p]=np.nanstd(Hallsims[:,(r_vals==1.1),skip:])
            
            #for i in range(3):
                #plt.plot(r_vals,Hallsims[i,:,19],color='black',linewidth=2, alpha=0.3)
                
        else:
            print("file " + file + " not found")
            
        #Get molecule number averages
        with open(folder+'.pkl', 'rb') as handle:
            data = pickle.load(handle)
        #data = { 'params':parameters,
#        'mtot42t': meantot42t, 'mtot42t2': meantot42t2,
#        'mgef42t': meangef42t, 'mgef42t2': meangef42t2,
#        'Nts':Ntimes, 'skp':skipframes, 'Ns':Nsims, 'smpl':sampling_sim,
#        }
        meanmolec[f,p]=data['mtot42t'][p]
        stddevmolec[f,p]= (np.asarray(data['mtot42t2'][p]) - np.asarray(data['mtot42t'][p])**2)**0.5
        stderrmolec[f,p]= stddevmolec[f,p]/((data['Nts'])*data['Ns'])**0.5
        
        
plt.close('all')

colors = [cm.hot(i) for i in np.linspace(0.6, 0, len(folders))]
matplotlib.rcParams.update({'font.size': 8})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1

fig=plt.figure(1)
        
for f in range(len(folders)):  
    #if f==0:
    #    plt.errorbar(params[1:],meanHs[f,1:],stdHs[f,1:],color=colors[f],alpha=1,marker="o",markersize=2,capsize=2,label=labels[f])       
    #else:
    plt.errorbar(params,meanHs[f,:],stdHs[f,:],color=colors[f],alpha=1,marker="o",markersize=2,capsize=2,label=labels[f])       
    
plt.title(title)
plt.ylabel(r'$H(r=1.1 \mu m)$')
plt.xlim([0,510])
plt.xlabel(r'GEF')
plt.xticks([0,100,200,300,400,500])
plt.legend(loc=2,fontsize=6)
fig.set_size_inches(2.5,2.3)
plt.ylim([-2,6])
#plt.axhline(y=1./60,color='black',linestyle='--')
#ax.set_xscale('log')
#ax.set_yscale('log')
plt.tight_layout()
fig.savefig(figurenameH)
plt.show()

#FIGURE OF MASS VS GEF
fig=plt.figure(2)

for f in range(len(folders)):  
    plt.errorbar(params,meanmolec[f,:],stddevmolec[f,:],color=colors[f],alpha=1,marker="o",markersize=2,capsize=2,label=labels[f])       

plt.title(title)
plt.ylabel(r'Total active Cdc42')
#plt.xlim([0,510])
plt.xlabel(r'GEF')
plt.xticks([0,100,200,300,400,500])
#plt.xlim([0,200])
#plt.ylim([0,500])
plt.legend(loc=2,fontsize=6)
fig.set_size_inches(2.5,2.3)
#plt.ylim([0,6])
#plt.axhline(y=1./60,color='black',linestyle='--')
#ax.set_xscale('log')
#ax.set_yscale('log')
plt.tight_layout()
fig.savefig(figurenameMass)
plt.show()

