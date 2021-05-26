function assemble_and_check_files(sourcedir,bashmetadatafile,outputfile)
% Assembles all the files from a single directory into a single file.
% Inputs:
%   sourcedir : the directory containing the files of interest. OK if other files are present in the same directory.
%   bashmetadatafile : the filename of the bashmetadata generated when running generate_bashfile().
%   outputfile : where to save the assembled results.

% bashmetadata
bmd = load(bashmetadatafile);

parallelChunkNames = cell(numel(1:bmd.chunksize:bmd.realizations),1);
results.t = cell(1,bmd.realizations);
results.A = cell(1,bmd.realizations);

expected_datapts = floor(bmd.nsteps/bmd.datagrain);

filecounter = 1;
for i=1:bmd.chunksize:bmd.realizations
    parallelChunkNames{filecounter} = [sourcedir bmd.matfile_prefix ...
                            strrep(sprintf('At%04i_Bt%04i_km%g_kr%g_DA%g_DB%g_realiz%04i_to_%04i',...
                                           bmd.Atot,bmd.Btot,bmd.k_micro,bmd.kr,bmd.DA,bmd.DB,i,i+bmd.chunksize-1),'.','-')];
    if exist([parallelChunkNames{filecounter},'.mat'])
        currdat = load(parallelChunkNames{filecounter});
        curridxs = i:(i+bmd.chunksize-1);
        results.t(curridxs(1):curridxs(end)) = currdat.t_TS;
        results.A(curridxs(1):curridxs(end)) = currdat.A_TS;
        invalidID = find(cellfun(@isempty,currdat.t_TS));
        if ~isempty(invalidID)
            fprintf('at least one entry was empty in %s\n\tentry:',parallelChunkNames{filecounter})
            disp(invalidID);
            for j=1:numel(invalidID)
                results.t{curridxs(invalidID(j))} = nan(expected_datapts,1);
                results.A{curridxs(invalidID(j))} = nan(expected_datapts,1);
            end
        end
    else
        curridxs = i:(i+bmd.chunksize-1);
        for j=curridxs(1):curridxs(end)
            results.t{j} = nan(expected_datapts,1);
            results.A{j} = nan(expected_datapts,1);
        end
    end
    filecounter = filecounter + 1;
end
save(outputfile,'results')
end
