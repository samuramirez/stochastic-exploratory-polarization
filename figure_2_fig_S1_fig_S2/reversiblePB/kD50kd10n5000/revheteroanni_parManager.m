function revheteroanni_parManager(NA,NB,DA,DB,kf,kr,sigma,xwin,nsteps,dt,datagrain,nRealizToRun,RNGstart,tempDir,parChunkName,numcpus)


mkdir(tempDir);
pc=parcluster('local');
pc.JobStorageLocation = tempDir;
parpool(pc,numcpus);

outputfile = parChunkName; % handle this externally so i can document in the bashfile generator

t_TS = cell(nRealizToRun,1);
A_TS = cell(nRealizToRun,1);
runtime = zeros(nRealizToRun,1);


parfor paridx=1:nRealizToRun
   currRealiz = RNGstart + paridx;
                           
   [t_TS{paridx},A_TS{paridx},runtime(paridx)] = revheteroanni_parCompat(NA,NB,DA,DB,kf,kr,sigma,xwin,nsteps,dt,datagrain,currRealiz)
end

save(outputfile,'t_TS','A_TS','runtime','NA','NB','DA','DB','kf','kr','sigma','xwin','dt','datagrain','nsteps','nRealizToRun','RNGstart');

end
