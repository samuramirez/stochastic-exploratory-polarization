import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io
import os.path
from matplotlib import cm

plt.ion()
matplotlib.rcParams.update({'font.size': 11})
#matplotlib.rcParams.update({'linewidth': 4})
matplotlib.rcParams['lines.linewidth'] = 1.5
#home = os.path.expanduser('~')

#kmicro over Deff

kaDPB = ["0p05", "0p5", "5","10", "50","100"]
kaDSG = ["0.05","0.50","5","10","50","100"]
#micro dissosiation coefficient
kdPB = ['0','0p1','1','10','100']
kdSG = ['0.0','0.1','1','10','100']

#kD 10

kaDind=4;
kdind=3;
nprot=5000
protr=0.005
nsims=1000 #10000 for nprot 5
Ns=[20,40,60]
L=1.0

number = len(Ns)
#cmap = plt.get_cmap('jet')
#colors = [cmap(i) for i in np.linspace(0, 1, number)]
#colors = [cm.jet(i) for i in np.linspace(0.8, 0.2, number)]
colors = [ cm.jet(0.2), cm.jet(0.5) ,cm.jet(0.9) ]
#linestyles=[(0,(1,0.5)),'dashed','solid']
linestyles=[(0,(1,0.5)),'solid','solid']

timeTS=0


#PLOT


plt.close('all')
#f, axarr = plt.subplots(2, 2)

#timeseries prot    
#matdata['results'][0][0][1][0][0]
#array of timeseries
#matdata['results'][0][0][1][0]
#rate = dA/(dt*meanA*meanB)
#meanA=(float(At)+mat['A_final'])/2

mets = ['ka','kh','ksrhet']
#mets = ['khyb']
metstitle=[r'$k_a$',r'$k_h$',r'$k_c$']
#metstitle=[r'$k_{hyb}$']

#mets = ['ksrhet']


#for figrow,kaDind in enumerate([2,3]):
    
#for figcol,kdind in enumerate([2,3]):

#LOAD PARTICLE BASED DATA
sim_dir = './reversiblePB/kD' + kaDPB[kaDind] + 'kd' + kdPB[kdind] + 'n' + str(nprot)
file =  sim_dir + '/assembled_data.mat'
if os.path.isfile(file):
    matdata = scipy.io.loadmat(file)
    #timeseries prot    
    #matdata['results'][0][0][1][0][0]
    #array of timeseries
    #matdata['results'][0][0][1][0]
    timeTSs=matdata['results'][0][0][0][0]
    timeTS=np.squeeze(timeTSs[0])
    Ats=matdata['results'][0][0][1][0]
    listAts = [np.squeeze(Ats[i]) for i in range(len(Ats))]
    Ats=np.asarray(listAts)
    Ats=np.squeeze(Ats)
    meanAt=np.mean(Ats,0)
    stdAtPB=np.std(Ats,0)
else:
    print("file " + file + " not found")
        
  
#figure index
fi=1    
    
for imet,met in enumerate(mets):  
    fig=plt.figure(fi)
    fi+=1
    
    #LOAD SPATIAL GILLESPIE DATA
    
    sim_dir ='./association5000molSG/'
    

    for iN, N in enumerate(Ns):
        filename= 'rev' + met + 'N' + str(N) + 'ka' + kaDSG[kaDind] + 'D' + 'kd' + kdSG[kdind] + 'nA0' + str(nprot)+'nB0'+str(nprot)+'.dat'
        file = sim_dir + filename 
        if os.path.isfile(file):
            #print(name)
            AtSG = np.genfromtxt(file,delimiter = ",")
            tSG = AtSG[:,0]
            meanAtSG = AtSG[:,1]
            mean2AtSG = AtSG[:,2]
            stdAtSG = (mean2AtSG - meanAtSG**2)**0.5
            #if iN==2:
             #   plt.fill_between(tSG,meanAtSG - stdAtSG, meanAtSG + stdAtSG, color=colors[iN], alpha=0.1)
                
            plt.plot(tSG,meanAtSG,label=r'$h =L/$'+ '{:.0f}'.format(N)+'='+'{:.1f}'.format(L/N/protr)+r'$\rho$'  , alpha=1,color=colors[iN],linestyle=linestyles[iN])
            #plt.errorbar(tSG,meanAtSG,stdAtSG,label='SG'+ met, alpha=0.6,color=cm.jet(1.0*iN/len(Ns)))
            plt.ylim(0,nprot)
            plt.xlim(0.0001,5)
            
            #axarr[figrow,figcol].set_xscale('log')

        else:
        	print("file " + file + " not found")
    #plt.errorbar(timeTS,meanAt,2*stdAt/nsims**0.5,label='PB', alpha=0.6, color='black')
    #plt.fill_between(timeTS,meanAt - stdAtPB, meanAt + stdAtPB, color='black', alpha=0.1)
    plt.plot(timeTS,meanAt,label='PB', alpha=1, color='black', linestyle='dashed')
    #plt.ylim(0,nprot)
    plt.xscale('log')

    plt.title(metstitle[imet])
    plt.xlabel('Time (s)')
    #plt.ylabel('# A')
    fig.set_size_inches(3.2, 2.7)
    #fig.set_size_inches(4, 3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    #fig.savefig('Ns_'+ metstitle[imet] +'_Np' + str(nprot) +'_kD' + kaDSG[kaDind] +'_kd_'+kdSG[kdind]+'.png', dpi=200)         


