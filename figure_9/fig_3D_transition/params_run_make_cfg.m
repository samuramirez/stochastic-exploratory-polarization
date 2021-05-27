filebase = '3d_newreac_gef100_k4a_1';
%0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 1
tstop = 1800;
dtime=1e-5;
samplingrate = 60*1/dtime;%60s
n_BemGEF = 100;
n_Cdc42 = 5000;
n_realiz = 30;
P4a = 1000*0.795774*dtime;

% The 'P's are bimolecular reaction probabilities, P=lambda_3d*dt.
% All the 'k's are first order rate constants, 1/s.
% The cytosol-to-membrane translocation rates (k1a, k5a) are scaled within
% make_smoldyn_cfg to the adsorption coefficients (kappa1a, kappa5a)  
% needed within Smoldyn, based on the domain's volume-to-surface-area ratio.

k1a = 0.1;
k1b = 10;
k2b = 0.63;
P2a = 25.46*dtime;
P3  = 55.70*dtime;
k4b = 10;
k5a = 4;
k5b = 6.5;
P7 = 8979.35*dtime;
P8 = 22448.39*dtime;
d_sphere = 4.5135;
rho = 0.02;
Dm = 0.0045;
Dc = 10;

% For setting up on Longleaf
%smoldyn_executable_path = '~/nas/longleaf/home/mikepab/smoldyn-2.56/cmake/smoldyn';
smoldyn_executable_path = '~/cmake/bin/smoldyn';
sh_name = sprintf('run_%s.sh',filebase);
fid=fopen(sh_name,'w'); % Create a bashfile to run the .cfg files.
fprintf(fid,'#!/bin/bash\n\n');

for i=1:n_realiz
  fileprefix = sprintf('%s_%02d',filebase,i);
  random_seed = i;

  % Creating configuration files
  make_smoldyn_cfg(fileprefix,tstop,dtime,samplingrate,random_seed,...
                   k1a,k1b,P2a,k2b,P3,P4a,k4b,k5a,k5b,P7,P8,...
                   Dm,Dc,rho,d_sphere,n_BemGEF,n_Cdc42);

  % Setting up the .sh file
  cfg_name = sprintf('%s.cfg',fileprefix);
  fprintf(fid,'sbatch -p general -N 1 -t 240:00:00 --wrap="%s %s"\n',smoldyn_executable_path,cfg_name);
end
fclose(fid);
