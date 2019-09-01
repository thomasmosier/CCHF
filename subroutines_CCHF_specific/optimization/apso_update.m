function [cfN, vel, lrnRt, state] = apso_update(cf, fitScore, prmBnds, lrnRt, statePrev, vel)
%Adaptive particle swarm optimization
%See description in "Zhi-Hui Zhan et al. - 2009 - Adaptive Particle Swarm
%Optimization.pdf"
%Full citation: Zhi-Hui Zhan, Jun Zhang, Yun Li, & Chung, H. S.-H. (2009). Adaptive 
%Particle Swarm Optimization. IEEE Transactions on Systems, Man, and 
%Cybernetics, Part B (Cybernetics), 39(6), 1362?1381. 
%https://doi.org/10.1109/TSMCB.2009.2015956


%Size of cf = [iterations, population, parameter]
%Size of fitScore = [iterations, population]

%Current solution for each member of population:
cfCurr = squeeze(cf(end,:,:));

%Size of current coefficients:
szCf = size(cfCurr);
if szCf(2) == 1 && szCf(1) ~= 1
    cfCurr = cfCurr';
    szCf = fliplr(szCf);
end

genCurr = numel(fitScore(:,1));

%Identify "personal best" (PBest) for each individual in the population 
[~, indPBest] = min(fitScore,[],1);
cfPBest = nan(szCf);
for ii = 1 : szCf(1)
    cfPBest(ii,:) = squeeze(cf(indPBest(ii),ii,:))';
end

%Find "global best" (GBest) for entire search algorithm:
if genCurr == 1
    genGBest = 1;
    [~, indGBest] = min(fitScore);
else
    [~, genGBest, indGBest] = min2d(fitScore);
end

if isempty(vel) || ~isequal(size(vel), szCf)
    warning('psoUpdate:velNotInitialized', ['The velocity matrix is ' ...
        'being set to zero because it is empty or not the correct size.']);
    vel = zeros(szCf); 
end

%Calculate distance between members (Eq. 7 in Zhi-Hui Zhan et a.):
d = zeros(szCf(1), 1);
indAll = (1:szCf(1));
for ii = 1 : szCf(1)
    indOther = setdiff(indAll,ii);
    for jj = 1 : szCf(1) - 1
        d(ii) = d(ii) + sqrt(nansum((cfCurr(ii,:) - cfCurr(indOther(jj),:)).^2));
    end
end
d = d/(sum(~isnan(cfCurr(:,1))) - 1);

%Calculate f (Eq. 8 in Zhi-Hui Zhan et a.)
[~, indDg] = min(fitScore(end,:));
dg = d(indDg);
dMx = nanmax(d);
dMn = nanmin(d);
f = (dg - dMn)/(dMx - dMn);

if f < 0 || f > 1
    error('psoUpdate:fOutOfRange', ['f is ' num2str(f) ', but the bounds should be [0, 1].'])
end

%Define state values:
%Exploration (Eq. 9a):
if f >= 0 && f <= 0.4
    ms1 = 0;
elseif f > 0.4 && f <= 0.6
    ms1 = 5*f - 2;
elseif f > 0.6 && f <= 0.7
    ms1 = 1;
elseif f > 0.7 && f <= 0.8
    ms1 = -10*f + 8;
elseif f > 0.8 && f <= 1
    ms1 = 0;
end

%Exploitation (Eq. 9b):
if f >= 0 && f <= 0.2
    ms2 = 0;
elseif f > 0.2 && f <= 0.3
    ms2 = 10*f - 2;
elseif f > 0.3 && f <= 0.4
    ms2 = 1;
elseif f > 0.4 && f <= 0.6
    ms2 = -5*f + 3;
elseif f > 0.6 && f <= 1
    ms2 = 0;
end

%Convergence (Eq. 9c):
if f >= 0 && f <= 0.1
    ms3 = 1;
elseif f > 0.1 && f <= 0.3
    ms3 = -5*f + 1.5;
elseif f > 0.3 && f <= 1
    ms3 = 0;
end

%Jumping Out (Eq. 9d):
if f >= 0 && f <= 0.7
    ms4 = 0;
elseif f > 0.7 && f <= 0.9
    ms4 = 5*f - 3.5;
elseif f > 0.9 && f <= 1
    ms4 = 1;
end

%Fuzzy logic (see Fig. 4)
if ms3 ~= 0 && ms2 == 0
    state = 3;
elseif ms3 ~=0 && ms2 ~= 0
    if statePrev == 1 || statePrev == 2 || ms2 > ms3
        state = 2;
    else
        state = 3;
    end
elseif ms2 ~= 0 && ms1 == 0
    state = 2;
elseif ms2 ~=0 && ms1 ~= 0
    if statePrev == 4 || statePrev == 1 || ms1 > ms2
        state = 1;
    else
        state = 2;
    end
elseif ms1 ~= 0 && ms4 == 0
    state = 1;
elseif ms1 ~=0 && ms4 ~= 0
    if statePrev == 3 || statePrev == 4 || ms4 > ms1
        state = 4;
    else
        state = 1;
    end 
elseif ms4 ~= 0
    state = 4;
else
    error('apsoUpdate:unknownLogicState', 'This logic state not expected.')
end

%Assign changes to the acceleration coefficients (c1 and c2) based on state determined above:
if genCurr == 1
    lrnRt = [2; 2];
else
    delta = randi([5,10])/100; %Paper recommends range of 0.05 to 0.1
    switch state 
        case 1
            lrnRt(1) = lrnRt(1) + delta;
            lrnRt(2) = lrnRt(2) - delta;
        case 2
            lrnRt(1) = lrnRt(1) + delta/2;
            lrnRt(2) = lrnRt(2) - delta/2;
        case 3
            lrnRt(1) = lrnRt(1) + delta/2;
            lrnRt(2) = lrnRt(2) + delta/2;
        case 4
            lrnRt(1) = lrnRt(1) - delta;
            lrnRt(2) = lrnRt(2) + delta;
        otherwise
            error('apsoUpdate:noLogicState', 'No logic state chosen.')
    end
    
    %Enforce min and max acceleration values:
    lrnMn = 1;
    lrnMx = 2.5;
    lrnRt(lrnRt < lrnMn) = lrnMn;
    lrnRt(lrnRt > lrnMx) = lrnMx;
    %Ensure that max sum is 4:
    if sum(lrnRt) > 4
        lrnRt = lrnRt * (4/sum(lrnRt));
    end
end

%Set the velocity inertia
if genCurr == 1
    omega = 0.9;
else
    omega = 1/(1 + 1.5*exp(-2.6*f));
end

if omega < 0.4
    omega = 0.4;
elseif omega > 0.9
    omega = 0.9;
end


%Update velocity:
vel = omega*vel ...
    + lrnRt(1) * rand(szCf) .* (cfPBest - cfCurr) ...
    + lrnRt(2) * rand(szCf) .* (repmat(squeeze(cf(genGBest, indGBest, :))', [szCf(1), 1]) - cfCurr);

%Limit vel to 20% search range:
vMax = 0.2;

sV = sign(vel);
rng = repmat((prmBnds(:,2) - prmBnds(:,1))', [szCf(1), 1]);
vel = sV.*min(abs(vel), vMax*rng);

%Create new parameters:
cfN = cfCurr + vel;

%If state 3, use elitist learning for the best performing parameter set:
if state == 3
    indPrm = randi([1,szCf(2)]);
    
    sdMx = 1;
    sdMn = 0.1;
    sdUse = sdMx - (sdMx - sdMn)*min(genCurr, 100)/100;
    
    [~, indEg] = max(fitScore(end,:));
    
    %Always replace the worst performer (this is slightly different than
    %reference)
    cfN(indEg, indPrm) = cfN(indEg, indPrm) + (prmBnds(indPrm,2) - prmBnds(indPrm,1))*randn(1)*sdUse;
end

%Ensure new coefficients are within parameter bounds:
for ii = 1 : szCf(2)
   indNeg = find(cfN(:,ii) < prmBnds(ii,1));
   indPos = find(cfN(:,ii) > prmBnds(ii,2));
   if ~isempty(indNeg)
       cfN(indNeg,ii) = prmBnds(ii,1);
   end
   if ~isempty(indPos)
       cfN(indPos,ii) = prmBnds(ii,2);
   end
end
    
    

%%FROM http://www.swarmintelligence.org/tutorials.php
% 3. The algorithm
% As stated before, PSO simulates the behaviors of bird flocking. Suppose 
%the following scenario: a group of birds are randomly searching food in an 
%area. There is only one piece of food in the area being searched. All the 
%birds do not know where the food is. But they know how far the food is in 
%each endation. So what's the best strategy to find the food? The 
%effective one is to follow the bird which is nearest to the food. 

% PSO learned from the scenario and used it to solve the optimization 
%problems. In PSO, each single solution is a "bird" in the search space. We 
%call it "particle". All of particles have fitness values which are 
%evaluated by the fitness function to be optimized, and have velocities 
%which direct the flying of the particles. The particles fly through the 
%problem space by following the current optimum particles. 
% 
% PSO is initialized with a group of random particles (solutions) and then 
%searches for optima by updating generations. In every endation, each 
%particle is updated by following two "best" values. The first one is the 
%best solution (fitness) it has achieved so far. (The fitness value is also 
%stored.) This value is called pbest. Another "best" value that is tracked 
%by the particle swarm optimizer is the best value, obtained so far by any 
%particle in the population. This best value is a global best and called 
%gbest. When a particle takes part of the population as its topological 
%neighbors, the best value is a local best and is called lbest.

% After finding the two best values, the particle updates its velocity and 
%positions with following equations:
% 
% v[] = v[] + c1 * rand() * (pbest[] - present[]) + c2 * rand() * (gbest[] - present[])
% present[] = persent[] + v[]
% 
% v[] is the particle velocity, persent[] is the current particle 
%(solution). pbest[] and gbest[] are defined as stated before. rand () is a 
%random number between (0,1). c1, c2 are learning factors. usually c1 = c2 = 2. 
% 
% The pseudo code of the procedure is as follows
% 
% For each particle 
%     Initialize particle
% END
% 
% Do
%     For each particle 
%         Calculate fitness value
%         If the fitness value is better than the best fitness value (pBest) in history
%             set current value as the new pBest
%     End
% 
%     Choose the particle with the best fitness value of all the particles as the gBest
%     For each particle 
%         Calculate particle velocity according equation (a)
%         Update particle position according equation (b)
%     End 
% While maximum endations or minimum error crendia is not attained
% 
% Particles' velocities on each dimension are clamped to a maximum velocity 
%Vmax. If the sum of accelerations would cause the velocity on that 
%dimension to exceed Vmax, which is a parameter specified by the user. Then 
%the velocity on that dimension is limited to Vmax.
