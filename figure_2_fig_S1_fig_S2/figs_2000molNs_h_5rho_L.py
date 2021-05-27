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

Lname=['L_0p4','L_0p8','L_1p6','L_3p2','L_6p4']
L = [0.4,0.8,1.6,3.2,6.4]
Lind=0

nprot=2000
protr=0.005
nsims=100 
Ns=[16,32,64,128,256]

N=Ns[Lind]

number = len(Ns)

colors =['magenta','cyan']
#linestyles=[(0,(1,0.5)),'dashed','solid']
linestyles=[(0,(1,0.5)),'solid','solid']

timeTS=0

alpha = float(kaDPB[kaDind])/(2*np.pi)
n=2000/np.asarray(Ns)**2

n=np.asarray([8,2,1,1,1])

Rc= 5*protr /(np.pi*n)**0.5
lamb = protr/Rc
#lamb=0.00000001
F= np.log(1/lamb)/(1-lamb**2)**2 - (3-lamb**2)/4/(1-lamb**2)
F*alpha

#Ns=[8,16,32,64,128]
#n=2000/np.asarray(Ns)**2

#n=np.asarray([31,8,2,1,1])

#Rc= 10*protr /(np.pi*n)**0.5
#lamb = protr/Rc
#F= np.log(1/lamb)/(1-lamb**2)**2 - (3-lamb**2)/4/(1-lamb**2)
#F*alpha


#PLOT


plt.close('all')


#mets = ['ka','kh','ksrhet']
mets = ['kh','ksrhet']

#metstitle=[r'$k_a$',r'$k_h$',r'$k_c$']
metstitle=[r'$k_h$',r'$k_c$']



#LOAD PARTICLE BASED DATA
sim_dir = './reversiblePB/kD' + kaDPB[kaDind] + 'kd' + kdPB[kdind] + 'n' + str(nprot) + '_' + Lname[Lind]
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
    #fi+=1
    
    #LOAD SPATIAL GILLESPIE DATA
    
    sim_dir ='./associationSG_h_5sigma_L/'    

    
    filename= 'rev' + met + 'N' + str(N) + 'ka' + kaDSG[kaDind] + 'D' + 'kd' + kdSG[kdind] + 'nA0' + str(nprot)+'nB0'+str(nprot)+'.dat'
    file = sim_dir + filename 
    if os.path.isfile(file):
        #print(name)
        AtSG = np.genfromtxt(file,delimiter = ",")
        tSG = AtSG[:,0]
        meanAtSG = AtSG[:,1]
        mean2AtSG = AtSG[:,2]
        stdAtSG = (mean2AtSG - meanAtSG**2)**0.5
      
        plt.plot(tSG,meanAtSG,alpha=1,color=colors[imet], label=metstitle[imet])
        #plt.errorbar(tSG,meanAtSG,2*stdAtSG,alpha=1,color=colors[imet], label=metstitle[imet])

        #plt.plot(tSG,meanAtSG,label=r'$h =L/$'+ '{:.0f}'.format(N)+'='+'{:.1f}'.format(L/N/protr)+r'$\rho$'  , alpha=1,color=colors[iN],linestyle=linestyles[iN])
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

plt.title("L="+str(L[Lind])+r'$\mu m$'+ r', h=5$\rho$'+', N='+str(N))
plt.xlabel('Time (s)')
#plt.ylabel('# A')
fig.set_size_inches(3.2, 2.7)
#fig.set_size_inches(4, 3)
plt.legend(fontsize=8)
plt.tight_layout()

#fig.savefig('N'+str(N) + '_h_5rho'+ '_Np' + str(nprot) +'_kD' + kaDPB[kaDind] +'_kd_'+kdPB[kdind]+'.pdf')         




