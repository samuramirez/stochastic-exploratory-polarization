DA=0.0025;
DB=0.0025;
ka = 50*2*DA;
kd = 0.0;
Atot = 5;
Btot= 5;
realizations=10000;
L=1;
sigma=0.005;
nsteps=10*10^5; 
dt=10^-5;
datagrain=10^2;
chunksize=1000;
matfile_prefix='rev';
bashfilename=['bash_' matfile_prefix];

generate_bashfile(DA,DB,ka,kd,Atot,Btot,realizations,L,sigma,nsteps,dt,datagrain,matfile_prefix,bashfilename,chunksize);