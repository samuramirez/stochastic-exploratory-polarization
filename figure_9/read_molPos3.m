function [t,positions] = read_molPos3(filename,nframes)

positions.Cdc42T = cell(nframes,1);
positions.BemGEF42 = cell(nframes,1);
positions.Cdc42Dm = cell(nframes,1);
positions.Cdc42Dc = cell(nframes,1);
positions.BemGEFm = cell(nframes,1);
positions.BemGEFc = cell(nframes,1);
positions.complex_Cdc42Dm_BemGEFm = cell(nframes,1);
positions.complex_Cdc42Dm_BemGEF42 = cell(nframes,1);
t=nan(nframes,1);
fid=fopen(filename);


frameid = 1;
while ~feof(fid)
    % each line contains:
    % t x1 y1 x2 y2 .. xn yn
    % for all n species of the particular time.
    % all 8 species are listed in the order
    % Cdc42T BemGEF42 Cdc42Dm Cdc42Dc BemGEFm BemGEFc complex_Cdc42Dm_BemGEFm complex_Cdc42Dm_BemGEF42
       
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    %currt
    %frameid
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.Cdc42T{frameid} = [x,y,z];
    
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.BemGEF42{frameid} = [x,y,z];
        
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);    
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.Cdc42Dm{frameid} = [x,y,z];
    
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.Cdc42Dc{frameid} = [x,y,z];
    
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.BemGEFm{frameid} = [x,y,z];
    
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.BemGEFc{frameid} = [x,y,z];
    
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    positions.complex_Cdc42Dm_BemGEFm{frameid} = [x,y,z];
    
    currline = fgetl(fid);
    [currt,x,y,z]=entry_to_xyz(currline);
    positions.complex_Cdc42Dm_BemGEF42{frameid} = [x,y,z];
    if ~isnan(currt) && isnan(t(frameid))
        t(frameid)=currt;
    end
    frameid=frameid+1;
end
fclose(fid);

end

function [t,x,y,z]=entry_to_xyz(line)
    coords=sscanf(line,'%g'); 
    if numel(coords)>1
        t = coords(1);
        x = coords(2:3:end);% skip the time value, which is the first entry
        y = coords(3:3:end);
        z = coords(4:3:end);
    else
        t = nan; x=[]; y=[]; z=[];
    end
end