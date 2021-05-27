function [t,steps]=calculate_activeCOM_steps(filename,nframes,R)
% calculate_activeCOM_steps: Generates a measure of patch wandering, "steps".
%   The center of mass for the active Cdc42 at each time point is calculated,
%   and the center of mass vector (starting from the origin) is projected to
%   the surface of the simulated sphere. The steps are the great circle
%   distances between individual projected COMs. A highly mobile patch will
%   take large steps. The time step chosen is an important factor in steps;
%   here we calculate between each recorded step (typically 10 sec).
%
%   Inputs:
%     filename: .xyz file with molecular coordinates
%     nframes: maximum # of frames that is in the file.
%     R: radius of the simulated sphere.
%   Outputs:
%     t: time points for steps.
%     steps: great circle distances between projected active Cdc42 COMs, as above.

% Get coordinates
[t,pos]=read_molPos3(filename,nframes);

% Get active Cdc42 center of mass coordinates
active_COMs = get_active_COMs(pos,nframes,R);

% Project the center of mass coordinates to the sphere surface
active_COMs_proj = project_COM_to_sphere_surface(active_COMs,R);

% Calculate great circle distances between vectors
steps=get_stepsizes(active_COMs_proj,R);

end

function [active_COMs]=get_active_COMs(positions,nframes,R)
% Returns ORIGIN-CENTERED (0,0,0) center-of-mass
% coordinates of active Cdc42 molecules

active_COMs = zeros(nframes,3);
for i=1:nframes
if ~isempty(positions.Cdc42T{i}) && ~isempty(positions.BemGEF42{i})
    all_active_x = [positions.Cdc42T{i}(:,1);positions.BemGEF42{i}(:,1)];
    all_active_y = [positions.Cdc42T{i}(:,2);positions.BemGEF42{i}(:,2)];
    all_active_z = [positions.Cdc42T{i}(:,3);positions.BemGEF42{i}(:,3)];
elseif isempty(positions.Cdc42T{i}) && ~isempty(positions.BemGEF42{i})
    all_active_x = positions.BemGEF42{i}(:,1);
    all_active_y = positions.BemGEF42{i}(:,2);
    all_active_z = positions.BemGEF42{i}(:,3);
elseif ~isempty(positions.Cdc42T{i}) && isempty(positions.BemGEF42{i})
    all_active_x = positions.Cdc42T{i}(:,1);
    all_active_y = positions.Cdc42T{i}(:,2);
    all_active_z = positions.Cdc42T{i}(:,3);
else
    all_active_x = nan;
    all_active_y = nan;
    all_active_z = nan;
end

tmp = [all_active_x-R,all_active_y-R,all_active_z-R];

active_COMs(i,:) = nanmean(tmp,1); % average along first dimension
end
end

function [adjusted_COMs]=project_COM_to_sphere_surface(COMs,R)
nframes=size(COMs,1);
adjusted_COMs=zeros(nframes,3);
for i=1:nframes
   scaleFac = R/norm(COMs(i,:));
   adjusted_COMs(i,:)=scaleFac*COMs(i,:);
end
end

function [d] = get_stepsizes(COMs,R)
nframes=size(COMs,1);
d=zeros(nframes,1);
for i=2:nframes
   [az1,el1,~]=cart2sph(COMs(i-1,1),COMs(i-1,2),COMs(i-1,3));
   [az2,el2,~]=cart2sph(COMs(i,1),COMs(i,2),COMs(i,3));

   d(i) = great_circle_distance(az1,el1,az2,el2,R);
end
end
