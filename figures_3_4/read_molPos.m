function positions = read_molPos(filename,nframes)

positions.Cdc42T = cell(nframes,1);
positions.BemGEF42 = cell(nframes,1);
positions.Cdc42Dm = cell(nframes,1);
positions.Cdc42Dc = cell(nframes,1);
positions.BemGEFm = cell(nframes,1);
positions.BemGEFc = cell(nframes,1);
positions.complex_Cdc42Dm_BemGEFm = cell(nframes,1);
positions.complex_Cdc42Dm_BemGEF42 = cell(nframes,1);

fid=fopen(filename);

frameid = 1;
while ~feof(fid)
    % each line contains:
    % t x1 y1 x2 y2 .. xn yn
    % for all n molecules of a particular species.
    % Species are listed as rows in in the order
    % Cdc42T BemGEF42 Cdc42Dm Cdc42Dc BemGEFm BemGEFc complex_Cdc42Dm_BemGEFm complex_Cdc42Dm_BemGEF42
    
    currline = fgetl(fid);
    [x,y]=entry_to_xy(currline);
    positions.Cdc42T{frameid} = [x,y];
    
    currline = fgetl(fid);
    [x,y]=entry_to_xy(currline);
    positions.BemGEF42{frameid} = [x,y];
        
    %currline = fgetl(fid);
    %[x,y]=entry_to_xy(currline);
    %positions.Cdc42Dm{frameid} = [x,y];
    
    %currline = fgetl(fid);
    %[x,y]=entry_to_xy(currline);
    %positions.Cdc42Dc{frameid} = [x,y];
    
    %currline = fgetl(fid);
    %[x,y]=entry_to_xy(currline);
    %positions.BemGEFm{frameid} = [x,y];
    
    %currline = fgetl(fid);
    %[x,y]=entry_to_xy(currline);
    %positions.BemGEFc{frameid} = [x,y];
    
    %currline = fgetl(fid);
    %[x,y]=entry_to_xy(currline);
    %positions.complex_Cdc42Dm_BemGEFm{frameid} = [x,y];
    
    %currline = fgetl(fid);
    %[x,y]=entry_to_xy(currline);
    %positions.complex_Cdc42Dm_BemGEF42{frameid} = [x,y];
    
    frameid=frameid+1;
end
fclose(fid);

end

function [x,y]=entry_to_xy(line)
    coords=sscanf(line,'%g'); 
    if numel(coords)>1      
        x = coords(2:2:end);% skip the time value, which is the first entry
        y = coords(3:2:end);
    else
        x=[]; y=[];
    end
end