function [H_vals,r_vals] = PB_Hfn_clustQuant(filename)
% ========================================================================
% Mike Pablo - Aug 23, 2016.
% Uses a rescaled version of a well-known clustering metric.
% The H-function is a modification of the established Ripley's K-function.
% It is defined using the cumulative distribution of inter-particle
% distances, 
%          / r
%  K(r) =  |    P(r')dr'
%          / 0
%
% where P(r)dr is the probability of finding a particle lying between r and
% r + dr.
%
% Then, H(r) is defined as
%           ___________
%          |  L^2
% H(r) =  \| ____ K(r)   - r
%          |  pi
%
% where the -r represents the contribution from a complete spatial random 
% distribution. Therefore, non-zero H(r) represents deviation from complete
% spatial randomness at the particular inter-particle distance r.
%
% See Wehrens M, Rein ten Wolde P, and Mugler A. J Chem Phys 2014, 141,
% 205102-1 for more details.
% ========================================================================


%% [BEGIN CLUSTER ANALYSIS] =============================================
disp('Beginning cluster analysis.')
% Parameters for PCF
%rho=0.02;
L=8;
xwin=L/2;
ywin=xwin;
%dr = rho/2;
dr=0.1;

%Import data
tseries1l=load(filename);%one vector for each timepoint 
N=sqrt(size(tseries1l,2));
Nf=size(tseries1l,1);
Ntime=Nf;
timepoints=[1:Ntime];
dx=L/N;

maxpart = 5000;


%put tseries1L into 3 dimensional matrix form, from data every 10 seconds (1, 10, 20 ...)
tseries =[];
for k=1:Nf
    tseries(k,:,:)=reshape(tseries1l(k,:), N, N);
end
		


%tseries=reshape(tseries1l(1,:), N,N);
%j=2;
%for k=10:10:Nf
%tseries(:,:,j)=reshape(tseries1l(k,:), N, N);
%j=j+1;
%end

timeanalysis=[1:Nf];


xst=-ones(maxpart,length(timeanalysis));
yst=-ones(maxpart,length(timeanalysis));

%Get positions of particles from grid based simulations
for k=1:length(timeanalysis)
	xs=[];
	ys=[];
	for i=1:N
		for j=1:N
			if tseries(k,i,j)>0
				xcell=(i-1)*dx + rand(1,tseries(k,i,j))*dx;
				ycell=(j-1)*dx + rand(1,tseries(k,i,j))*dx;	
				xs=[xs xcell];
				ys=[ys ycell];
			end
		end
    end
   % display(k)
   % display(xs)
	xst(1:length(xs),k)=xs;
	yst(1:length(ys),k)=ys;
end	

xst(xst==-1)=nan;
yst(yst==-1)=nan;					

%max_r  = sqrt( (xwin*2)^2 + (ywin*2)^2 )/2;
max_r  = 5;
r_vals = 0:dr:max_r; % Will stop such that max(max(r_vals)) <= max_r
							  
PCFs   = zeros(length(r_vals),length(timeanalysis));

H_vals = zeros(length(r_vals),length(timeanalysis));
if xwin ~= ywin
    error('Warning - definitions used for H-function are not valid if xwin =/= ywin');
end

for i=1:length(timeanalysis)
    [PCFs(:,i),~] = computePr(xst(:,i),yst(:,i),dr,xwin,ywin);
    for j=1:length(r_vals)
        if j==1
            H_vals(j,i) = 0;
        else
            H_vals(j,i) = sqrt( (2*xwin)^2 / pi * trapz(r_vals(1:j),PCFs(1:j,i))) - r_vals(j);
        end
    end
    fprintf('Computed %i of %i timesteps\n',i,length(timeanalysis));
end                                            

end %end function [H_vals,r_vals] = PB_Hfn_clustQuant(filename)


function [P,r_vals] = computePr(Xs,Ys,dr,xwin,ywin)
%disp('Computing a variation on the pair-wise distance distribution function')
% Compute a variation on the pair-wise distance distribution function P(r).
%          1    
% P(r) = -----  * sum_from_i=1_to_N {m_i(r)}
%        N(N-1)
% where m_i(r)dr is the number of particles found between r+dr.
% ===========================================================
% INPUTS
%   Xs   : x-coordinates of all particles of interest
%   Ys   : y-coordinates of all particles of interest
%   dr   : granularity of radial binning (same units as coordinates)
%   xwin : half-window length of domain (same units as coordinates)
%   ywin : half-window length of domain (same units as coordinates)
% OUTPUTS
%   P    : distribution of pairwise distances, also called the pair
%          correlation function. Units 1 / distance.
%
% This does take the periodic boundaries into consideration.
% Clusters spread across a boundary will be handled almost-correctly, in
% that they get double-counted if you go out too far.
% ===========================================================

% Xs and Ys MUST be in Nx2 format, where N is the # of particles
Xs = Xs(~isnan(Xs)); Ys = Ys(~isnan(Ys));
coordMat = [Xs Ys];
[r,c] = size(coordMat);
if c ~= 2
    disp('This function is meant for 2-dimensional data, use at own risk.')
end
Nparts = r;
sqDists = pdist(coordMat,@(ZI,ZJ) SQeuclTorus(ZI,ZJ,xwin*2,ywin*2));
sqDistsMat=squareform(sqDists);
% The rows of sqDists. Compute out to the half-diagonal of the domain.
%max_r  = sqrt( (xwin*2)^2 + (ywin*2)^2 )/2;
max_r  = 5;
r_vals = 0:dr:max_r; % Will stop such that max(max(r_vals)) <= max_r
r2_vals = dr:dr:(max_r+dr); % secondary r_vals to reduce parallel broadcast overhead

m_i = zeros(length(r_vals),1);
parfor r_idxs = 1:(length(r_vals)-1) % For each radius
    for i = 1:length(Xs) % For each particle
        for j = 1:length(Xs) % Check through all particles and see if they are close enough
            if i~=j
                
                %currsqDist = sqDists( (i-1)*(Nparts-i/2)+j-i ); %index
                %goes to 0 when Npart = 2. !
                currsqDist = sqDistsMat(i,j); 
                if currsqDist > r_vals(r_idxs)^2 && currsqDist < r2_vals(r_idxs)^2
                    m_i(r_idxs) = m_i(r_idxs) + 1; % Increase the count.
                end
            end 
        end
    end
end

% Normalize.
P = 1 / (Nparts * (Nparts - 1)) .* m_i / dr;

end

function [d2] = SQeuclTorus(ZI,ZJ,FULLXWIN,FULLYWIN)
% Calculates squared euclidean distances between points on a periodic 2D plane.
% The plane is described xwin, ywin.
% the inputs:
%     ZI - some arbitrary 1-by-n vector that contains a single observation
%     from X or Y in pdist2(X,Y,@euclTorus(ZI,ZJ) euclTorus(X,Y,XWIN,YWIN))
%     ZJ - the m2-by-n matrix containing multiple observations from Y or X.
%the output:
%     d2 - an m2-by-1 vector of distances whose Jth element is the distance
%     between the observations ZI and ZJ(J,:).

diff=bsxfun(@minus,ZI,ZJ); % subtract each coordinate in the ZI observation...
                     % from each coordinate in the ZJ(j)-th observation.
                     % Return a matrix size ZJ.

% Take the square room of the sum of the minimized squares
d2=bsxfun(@min,abs(diff(:,1)),FULLXWIN-abs(diff(:,1))).^2 + ...
   bsxfun(@min,abs(diff(:,2)),FULLYWIN-abs(diff(:,2))).^2;

end

% Function to establish data structures.
function [outfile,todays_date]=SetFolderbyDate(filebasename)
if ispc
    todays_date = num2str(yyyymmdd(datetime)); % folder
    todays_date = strcat(todays_date,'\');
    outfile     = strcat(todays_date,'\',filebasename);
elseif ismac
    todays_date = num2str(yyyymmdd(datetime)); % folder
    todays_date = strcat(todays_date,'/');
    outfile     = strcat(todays_date,'/',filebasename);
elseif isunix
    % Minor difference from actual dirs in the *nix format..
    % pc and mac give YYYYMMDD
    % *nix gives YYYY_M_D, where M and D can have 1 or 2 digits. Added
    % underscores for clarity; just note you cannot integrate
    % *nix names directly with pc or mac names if automating analysis
    mytime=clock;
    mytime=mytime(1:3); % get YYYYMD only
    todays_date=strrep(regexprep(num2str(fix(mytime)),' +',' '),' ','_');
    todays_date = strcat(todays_date,'/');
    outfile     = strcat(todays_date,'/',filebasename);
end
end


