import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io
import os.path
from matplotlib import cm

plt.ion()
matplotlib.rcParams.update({'font.size': 16})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1.5

Nsims=30;

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
Nmeth=["40","80","160"];

#FOR FIGURE 4B
meth=['kh','ksrhet'];
colors = ['black', cm.jet(0.3), 'red']
idxNmetplot=1
figurename='boxplotkhksrhet_h5rho.pdf'


#FOR FIGURE 4D
#meth=['ksrhet'];
#colors = ['black','red']
#idxNmetplot=2
#figurename='boxplotksrhet_h2p5rho.pdf'


#Nmeth=["120"];

GEFconc=["100","200","300","400","500","600"];
#GEFconc=["600","700"];
methodname=[r'$k_h$',r'$k_c$']


#colors = [cm.jet(i) for i in np.linspace(1, 0, len(GEFconc))]
#colNs= [cm.jet(i) for i in np.linspace(1, 0, len(Nmeth))] 
#cols= [ cm.jet(0.2) ,'red']
#markers=['x','^']
linestyles=['dotted','dashed','solid']

fi=1

  
# Clustering as a function of GEF averaging the last time points (250-300s) : H(r==1.1) 
#fig=plt.figure(fi)
#fi+=1
#PARTICLE BASED
#HvsGEFPB=[]
#HvsGEFPBstddev=[]
HvsGEFPBbox=np.zeros((len(GEFconc),Nsims))

for idxgef,gefconc in enumerate(GEFconc):
    file =  filename='./HvaluesPB/Hvaluesgef' + gefconc + 'PB.mat';
    if os.path.isfile(file):
        matdata = scipy.io.loadmat(file)        
        Hallsims=matdata['Hallsims']
        #HvsGEFPB.append(np.nanmean(Hallsims[:,(r_vals==1.1),-5:],(0,2)))
        #HvsGEFPBstddev.append(np.nanstd(Hallsims[:,(r_vals==1.1),-5:],(0,2)))
        #HvsGEFPBbox[idxgef,:]=Hallsims[:,(r_vals==1.1),-5:].flatten()
        
        for sim in range(Nsims):
            HvsGEFPBbox[idxgef,sim]= np.nanmean(Hallsims[sim,(r_vals==1.1),-5:])
    else:
        print("file " + file + " not found")        
#SPATIAL GILLESPIE
#HvsGEFSG=np.zeros((len(meth),len(Nmeth),len(GEFconc)))
#HvsGEFSGstddev=np.zeros((len(meth),len(Nmeth),len(GEFconc)))
HvsGEFSGbox=np.zeros((len(meth),len(Nmeth),len(GEFconc),Nsims))
for idxgef,gefconc in enumerate(GEFconc):       
    for idxmet,met in enumerate(meth):
        for idxNmet,Nmet in enumerate(Nmeth):
            file =  filename='./HvaluesSG/Hvalues'+'gef'+ gefconc +met+'N'+Nmet+'.mat';
            if os.path.isfile(file):
                matdata = scipy.io.loadmat(file)        
                Hallsims=matdata['Hallsims']      
                #HvsGEFSG[idxmet,idxNmet,idxgef]=np.nanmean(Hallsims[:,(r_vals==1.1),-5:],(0,2))
                #HvsGEFSGstddev[idxmet,idxNmet,idxgef]=np.nanstd(Hallsims[:,(r_vals==1.1),-5:],(0,2))
                #HvsGEFSGbox[idxmet,idxNmet,idxgef,:]=Hallsims[:,(r_vals==1.1),-5:].flatten()
                for sim in range(Nsims):
                    HvsGEFSGbox[idxmet,idxNmet,idxgef,sim]= np.nanmean(Hallsims[sim,(r_vals==1.1),-7:])    
            else:
                print("file " + file + " not found")
#MAKE PLOTS

fig, axes = plt.subplots(nrows=1, ncols=len(GEFconc), figsize=(15, 4))
bplots=[]
for idxgef,gefconc in enumerate(GEFconc):       
    #Hs=[HvsGEFPBbox[idxgef,:][(~np.isnan(HvsGEFPBbox[idxgef,:]))] ]
    Hs=[HvsGEFPBbox[idxgef,:][(HvsGEFPBbox[idxgef,:]>0)] ]
    for idxmet,met in enumerate(meth):
        #Hs.append(HvsGEFSGbox[idxmet,idxNmetplot,idxgef,:][(~np.isnan(HvsGEFSGbox[idxmet,idxNmetplot,idxgef,:]))])
        Hs.append(HvsGEFSGbox[idxmet,idxNmetplot,idxgef,:][(HvsGEFSGbox[idxmet,idxNmetplot,idxgef,:]>0)])
        #Hs.append(HvsGEFSGbox[idxmet,idxNmetplot,idxgef,:]) 
        
    bplot = axes[idxgef].boxplot(Hs,
                         vert=True,  # vertical box alignment
                         patch_artist=True,  # fill with color
                         #labels=labels,
                         widths=0.8,
                         showfliers=False
                         )
    bplots.append(bplot)
    for bplot in bplots:
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_color(color)
        for line, color in zip(bplot['medians'], colors):
            line.set_color('white')
        for idxline,line in enumerate(bplot['whiskers']):
            line.set_color(colors[int(idxline/2)])
        for idxline,line in enumerate(bplot['caps']):
            line.set_color(colors[int(idxline/2)])
            

for idxgef,ax in enumerate(axes):
    #ax.axis('off')
    #x positions of boxplots are integers
    ax.set_xlim(0,4)
    ax.set_ylim(0,3.5)

    #ax.yaxis.grid(True)
    #ax.set_xlabel('Three separate samples')
    ax.set_xlabel(GEFconc[idxgef])
    

# Only show ticks on the left and bottom spines
axes[0].axis('on')
axes[0].yaxis.set_ticks_position('left')
axes[0].set_xticks([])
axes[0].set_ylabel('clustering')

# Hide the right and top spines
axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)
#axes[0].spines['bottom'].set_visible(False)

#set position [left, bottom, width, height]
#get positions of first subplot
pos0=axes[0].get_position().bounds
#decrease width of first subplot by a factor
axes[0].set_position(pos0*np.array([1,1,0.4,1]))
#get new positions of first subplot
pos0=axes[0].get_position().bounds

#set positions of next subplots similar as the first but shifting to the right
for i in range(1,len(GEFconc)):
    axes[i].set_position(pos0+np.array([i*0.1,0,0,0]))
    axes[i].set_xticks([])
    axes[i].set_yticks([])
    axes[i].spines['right'].set_visible(False)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['left'].set_visible(False)

#axes[1].set_position([0.125, 0.125, 0.477272727273, 0.88])

#plt.xlabel(r'Total BemGEF')
#plt.ylabel(r'H(r=1.1$\mu m$,t=250-300s)')
#plt.ylabel(r'polarization')
#plt.ylim(0,3.5)
#plt.xlim(80,750)
fig.set_size_inches(6, 4)
#fig.set_size(5, 4)
#plt.legend(loc=4,fontsize=7.5)
plt.tight_layout()              
#fig.savefig(figurename)

'''
# Random test data
np.random.seed(19680801)
all_data = [np.random.normal(0, std, size=100) for std in range(1, 4)]
labels = ['x1', 'x2', 'x3']

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
bplots=[]
# rectangular box plot
bplot1= axes[0].boxplot(all_data,
                         vert=True,  # vertical box alignment
                         patch_artist=True,  # fill with color
                         labels=labels,
                         widths=0.5)  # will be used to label x-ticks
#axes[0].set_title('Rectangular box plot')

# notch shape box plot
bplot2 = axes[1].boxplot(all_data,
                         #notch=True,  # notch shape
                         vert=True,  # vertical box alignment
                         patch_artist=True,  # fill with color
                         labels=labels,
                         widths=0.5)  # will be used to label x-ticks
#axes[1].set_title('Notched box plot')

# fill with colors
colors = ['black', 'blue', 'red']
for bplot in (bplot1, bplot2):
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
#axes[0].set_xlim(0,8)


# adding horizontal grid lines
for ax in axes:
    ax.axis('off')
    ax.set_xlim(0,4)
    #ax.yaxis.grid(True)
    #ax.set_xlabel('Three separate samples')
    #ax.set_ylabel('Observed values')



# Only show ticks on the left and bottom spines
axes[0].axis('on')
axes[0].yaxis.set_ticks_position('left')
axes[0].set_xticks([])
# Hide the right and top spines
axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)
axes[0].spines['bottom'].set_visible(False)
#set position [left, bottom, width, height]
#get positions of first subplot
pos0=axes[0].get_position().bounds
#decrease width of first subplot by a factor
axes[0].set_position(pos0*np.array([1,1,0.2,1]))
#get new positions of first subplot
pos0=axes[0].get_position().bounds
#set positions of second subplot similar as the first but shifted to the right
axes[1].set_position(pos0+np.array([0.07,0,0,0]))
#axes[1].set_position([0.125, 0.125, 0.477272727273, 0.88])


plt.show()
'''