function [t_TS,A_TS,runtime] = revheteroanni_parCompat(NA,NB,DA,DB,kf,kr,sigma,xwin,nsteps,dt,datagrain,rngseed)
%Bimolecular simulation of A + A -> 0 with Yogurtcu et al theory
% ======================================================================

tic

rng(rngseed);

%% -------------------------
% Initialize pre-calculated values.
% Initialize data matrices.
%% --------------------------
% Compute lambda from 2D rate constant kf.
prxn=kf*dt/(sigma^2*pi);
% compute reverse rxn probability
prxn_rev = 1-exp(-kr*dt);

ywin=xwin;
As = zeros(NA,3); % NA x {x,y,state} % for the sake of the reversible rxns, let 0 be unbound, and integer values represent the bound particle index.
Bs = zeros(NB,3); % NB x {x,y,state}

As(:,1) = unifrnd(-xwin,xwin,NA,1);
As(:,2) = unifrnd(-ywin,ywin,NA,1);
As(:,3) = 0;

Bs(:,1) = unifrnd(-xwin,xwin,NB,1);
Bs(:,2) = unifrnd(-ywin,ywin,NB,1);
Bs(:,3) = 0;

A_TS = zeros(floor(nsteps/datagrain),1);
t_TS = zeros(floor(nsteps/datagrain),1);

%% ------------------------
% Begin iterative simulation
%% ------------------------
%i=1; % time index
for i=1:nsteps
    % Get annihilation reaction list and handle flags.
    [A_annihilated,B_annihilated] = getReactionList(As,Bs,sigma,prxn,xwin);
    
    % reverse reaction
    [A_dissoc,B_dissoc]=getFirstOrderRxnList(As,prxn_rev);
    
    % settle flags
    inertAs = true(NA,1);
    inertBs = true(NB,1);
    inertAs(A_annihilated) = false;
    inertAs(A_dissoc) = false;
    inertBs(B_annihilated) = false;
    inertBs(B_dissoc) = false;
    
    As(A_annihilated(:),3) = B_annihilated;
    Bs(B_annihilated(:),3) = A_annihilated;
    As(A_dissoc,3) = 0;
    Bs(B_dissoc,3) = 0;
        
    % Calculate diffusional updates
    freeA = ~logical(As(:,3)) & inertAs;
    As(freeA,1) = calcPos(As(freeA,1),DA,dt,xwin);
    As(freeA,2) = calcPos(As(freeA,2),DA,dt,ywin);
    
    freeB = ~logical(Bs(:,3)) & inertBs;
    Bs(freeB,1) = calcPos(Bs(freeB,1),DB,dt,xwin);
    Bs(freeB,2) = calcPos(Bs(freeB,2),DB,dt,ywin);
    
    % Displace dissociated As from their B
    dissocFlagA = false(NA,1);
    dissocFlagA(A_dissoc) = true;
    if any(dissocFlagA)
        [newAposx,newAposy] = CalcPos_Dissociation(dissocFlagA,As(:,1),As(:,2),xwin,ywin,sigma);
        As(:,1) = newAposx;
        As(:,2) = newAposy;
    end
    % Diffuse bound As, and move their bound B with them.
    boundA = logical(As(:,3));
    As(boundA,1) = calcPos(As(boundA,1),DA,dt,xwin);
    As(boundA,2) = calcPos(As(boundA,2),DA,dt,ywin);
    followerB = As(boundA,3);
    Bs(followerB,1) = As(boundA,1);
    Bs(followerB,2) = As(boundA,2);

    if mod(i-1,datagrain) == 0
        grainstep = ceil(i/datagrain);
        A_TS(grainstep) = sum(~boundA); % free A
        t_TS(grainstep) = i*dt;
    end
end


%% ------------------------------------------------
% Timestep iteration over.
%% ------------------------------------------------
runtime = toc;
end


% Function to find paired association reactions
function [listA,listB] = getReactionList(partsA,partsB,rad,prxn,xwin)
% Compute reacting pairs. If multiple potential reactants, checks in order of
% the particles' indices, since we only allow 1-1 stoichiometry. Returns flagged particle indices.

bool_reactive_A = ~logical(partsA(:,3)); % Require non-annhilated 
bool_reactive_B = ~logical(partsB(:,3)); % Require non annihilated


if sum(bool_reactive_A) > 0 && sum(bool_reactive_B) > 0 % More than one reactant needed
    posA=partsA(bool_reactive_A,1:2);
    posB=partsB(bool_reactive_B,1:2);
    % Calculate pairwise distances based on periodic 2D plane (torus)
    % Efficiency - may be better to call pdist with anonymous function.
    SQeuclDist = sqEuclTorus_mex(posA,posB,xwin*2);
    ReactionMat = SQeuclDist < rad^2;
    ReactionMat(ReactionMat) = rand(sum(sum(ReactionMat)),1) < prxn;

    [a1,a2] = find(ReactionMat);
    
    if ~isempty(a1) % if a1 is empty, a2 must be as well.
        listA_prefilt = zeros(length(a1),1);
        listB_prefilt = zeros(length(a2),1);
        translate_reactive_A = find(bool_reactive_A);
        translate_reactive_B = find(bool_reactive_B);
        
        listA_prefilt(:) = translate_reactive_A(a1);
        listB_prefilt(:) = translate_reactive_B(a2);

        % Check uniqueness of binding
        a_multiplicity = sum(ReactionMat,2); % Along rows
        b_multiplicity = sum(ReactionMat,1); % Along columns

        if any(a_multiplicity > 1) || any(b_multiplicity > 1)
            % Filter out duplicated indices
            [~,uniqueA] = unique(listA_prefilt);
            [~,uniqueB] = unique(listB_prefilt);
            oneToOneFilter = intersect(uniqueA,uniqueB);
            
            listA = listA_prefilt(oneToOneFilter);
            listB = listB_prefilt(oneToOneFilter);
        else % No need to filter, no multiplicities.
            listA = listA_prefilt;
            listB = listB_prefilt;
        end
    else
        listA = []; % If no reactives found after prxn check
        listB = [];
    end
else
    listA = []; % If no reactives found after the initial reactives check
    listB = [];
end
end

function [listA,listB] = getFirstOrderRxnList(partsA,revprob)
bool_reactive_A = logical(partsA(:,3)); % bound particles
doRxn = rand(sum(bool_reactive_A),1)<revprob;
candidatesA=find(bool_reactive_A);
% dissociate the molecules
listA = candidatesA(doRxn);
listB = partsA(listA,3);
end

% Function to perform a diffusion step
function [newPos] = calcPos(pos, D, dt, bound)
% Calculate new position
% The position function is updated (symetrically for x and y)
% as x_i(t+dt)=x_i(t)+sqrt(2D*dt)Z_i
% where
% x_i is an arbitrary coordinate (x1 = x, x2 = y)
% t is time
% dt is the timestep
% D is the diffusion coefficient of the particle
% Z_i is an arbitrary Gaussian random variable, mean=0, stdev=1
% ** pre-calculates the list of random variables at start of timestep
% for all particles for each step, then pass one here
% to minimize calls to normrnd
% note -- here, pos == x_i(t)
%               z   == Z_i

z=normrnd(0,1,size(pos));
delta=sqrt(2*D*dt)*z;
newPos=pos+delta;
newPos=adjustForBounds(newPos,bound); 

end

% Compute new positions for all particles of A after dissociation from
% complexes AB. The B particles do not move.
% Expects oldApos as a Nx1 matrix in 1D, or Nx2 matrix in 2D.
function [newAposx,newAposy] = CalcPos_Dissociation(dissocFlagA,oldAposx,oldAposy,xwin,ywin,sigma)
randTheta = rand(length(dissocFlagA),1)*2*pi;

% In 2D
newAposx = oldAposx + sigma*cos(randTheta).*logical(dissocFlagA);
newAposy = oldAposy + sigma*sin(randTheta).*logical(dissocFlagA);

newAposx = adjustForBounds(newAposx,xwin);
newAposy = adjustForBounds(newAposy,ywin);
end

% Function to adjust the input coordinates for the boundary conditions.
% Can use as oldXY == newXY to see if a change occurs
function [pos] = adjustForBounds(pos, bound)

% PERIODIC BOUNDARIES MODE FOR VECTOR pos
timesExceeded=floor((abs(pos)+bound)/(2*bound));
pos = pos -sign(pos)* 2*bound.*timesExceeded;

end
