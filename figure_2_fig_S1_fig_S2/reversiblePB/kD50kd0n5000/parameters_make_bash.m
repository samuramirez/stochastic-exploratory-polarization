DA=0.0025;
DB=0.0025;
ka = 50*2*DA;
kd = 0;
Atot = 5000;
Btot= 5000;
realizations=100;
L=1;
sigma=0.005;
nsteps=5*10^5; 
dt=10^-5;
datagrain=10^1;
chunksize=10;
matfile_prefix='rev';
bashfilename=['bash_' matfile_prefix];

generate_bashfile(DA,DB,ka,kd,Atot,Btot,realizations,L,sigma,nsteps,dt,datagrain,matfile_prefix,bashfilename,chunksize);
