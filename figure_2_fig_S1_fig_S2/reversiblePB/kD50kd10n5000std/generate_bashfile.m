function generate_bashfile(DA,DB,k_micro,kr,Atot,Btot,realizations,L,sigma,nsteps,dt,datagrain,matfile_prefix,bashfilename,chunksize)
% Generate a bash file to run homoannihilation reactions, and save the time series
% Inputs: 
%   DA,DB: diffusivity, um2/s
%   k_micro : microscopic reaction rate, um2/s
%   kr : reverse reaction rate, 1/s
%   Atot,Btot : total protein number
%   realizations : number replicates
%   L : full domain length (um); domain is LxL square
%   sigma : reaction radius (um)
%   dt : timestep (s)
%   datagrain : specifies how often data is recorded (=1, record each timestep; =10, record every 10 timesteps, etc.)
%   fractionCutoffA : fraction of A, below which the simulation can end.
%   matfile_prefix : the prefix for the outputted mat files
%   chunksize : the # of realizations to allot per node.

debug_mode = false;
if debug_mode
warning('Running in debug mode. Settings are probably not good for actual simulation.')
end

% input checks
assert(isequal(round(Atot),Atot) && Atot>0,'(Atot==%g) require Atot to be a non-negative integer',Atot)
assert(isequal(round(Btot),Btot) && Btot>0,'(Btot==%g) require Btot to be a non-negative integer',Btot)
assert(isequal(round(realizations),realizations) && realizations>0,'(realizations==%g) require realizations to be a non-negative integer',realizations)
assert(DA>0,'(DA==%g) require diffusivity to be positive',DA)
assert(DB>0,'(DB==%g) require diffusivity to be positive',DB)
assert(k_micro>=0,'(k_micro==%g) require reaction rate to be at least zero',k_micro)
assert(kr>=0,'(kr==%g) require reaction rate to be at least zero',kr)
assert(L>0,'(L==%g) require domain size to be positive',L)
assert(dt>0,'(dt==%g) require timestep to be positive',dt)
assert(sigma>0,'(sigma==%g) require sigma to be positive',sigma)
assert(isequal(round(datagrain),datagrain) && datagrain>0,'(datagrain==%g) require datagrain to be a non-negative integer',datagrain)
assert(isequal(round(chunksize),chunksize) && chunksize>0,'(chunksize==%g) require chunksize to be a non-negative integer',chunksize)
assert(isequal(round(nsteps),nsteps) && nsteps>0,'(nsteps==%g) require nsteps to be a non-negative integer',nsteps)


if (sqrt((DA+DB)*dt)/sigma)>0.1
    warning('the RMS step to sigma ratio is %g, perhaps shrink the timestep',sqrt((DA+DB)*dt)/sigma);
end

FID = fopen([bashfilename '.sh'],'w');
fprintf(FID,'#!/bin/bash\n\n');

tempParDirCounter = 1;

for realiz_idx = 1:chunksize:realizations
    tempParDirName = sprintf('tempDir%04i',tempParDirCounter);

    parallelChunkName = [matfile_prefix strrep(sprintf('At%04i_Bt%04i_km%g_kr%g_DA%g_DB%g_realiz%04i_to_%04i',...
                                                       Atot,Btot,k_micro,kr,DA,DB,realiz_idx,realiz_idx+chunksize-1),'.','-')];

    if debug_mode
        fprintf(FID,'sbatch -p general -N 1 -n 2 -t 00:10:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r ');
        fprintf(FID,'revheteroanni_parManager\\(%i,%i,%g,%g,%g,%g,%g,%g,%g,%i,%i,%i,%i,\\''%s\\'',\\''%s\\'',2\\) -logfile %04i.log"\n',...
                    Atot,Btot,DA,DB,k_micro,kr,sigma,L/2,nsteps,dt,datagrain,chunksize,realiz_idx-1,tempParDirName,parallelChunkName,tempParDirCounter);
    else
        fprintf(FID,'sbatch -p general -N 1 -n 12 -t 144:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r ');
        fprintf(FID,'revheteroanni_parManager\\(%i,%i,%g,%g,%g,%g,%g,%g,%g,%i,%i,%i,%i,\\''%s\\'',\\''%s\\'',12\\) -logfile %04i.log"\n',...
                    Atot,Btot,DA,DB,k_micro,kr,sigma,L/2,nsteps,dt,datagrain,chunksize,realiz_idx-1,tempParDirName,parallelChunkName,tempParDirCounter);
    end


    tempParDirCounter = tempParDirCounter + 1;
end
fclose(FID);
fprintf('Wrote to %s\n',bashfilename);
bash_metadata_filename = [bashfilename '_metadata'];
save(bash_metadata_filename,'DA','DB','k_micro','kr','Atot','Btot','realizations','L','sigma','dt','datagrain','nsteps','matfile_prefix','bashfilename','chunksize');
fprintf('Saved metadata at %s\n',bash_metadata_filename);

end
