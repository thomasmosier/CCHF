function [cfN, vel] = pso_update(cf, fitScore, prmBnds, lrnRt, vel)

%See description in "Zhi-Hui Zhan et al. - 2009 - Adaptive Particle Swarm
%Optimization.pdf"
%Full citation: Zhi-Hui Zhan, Jun Zhang, Yun Li, & Chung, H. S.-H. (2009). Adaptive 
%Particle Swarm Optimization. IEEE Transactions on Systems, Man, and 
%Cybernetics, Part B (Cybernetics), 39(6), 1362?1381. 
%https://doi.org/10.1109/TSMCB.2009.2015956


%and at bottom of this function. 


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

%Identify "personal best" (PBest) for each individual in the population 
[~, indPBest] = min(fitScore,[],1);
cfPBest = nan(szCf);
for ii = 1 : szCf(1)
    cfPBest(ii,:) = squeeze(cf(indPBest(ii),ii,:))';
end
%Find "global best" (GBest) for entire search algorithm:
%Find "global best" (GBest) for entire search algorithm:
if numel(fitScore(:,1)) == 1
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


%Update velocity:
vel = vel ...
    + lrnRt * rand(szCf) .* (cfPBest - cfCurr) ...
    + lrnRt * rand(szCf) .* (repmat(squeeze(cf(genGBest, indGBest, :))', [szCf(1), 1]) - cfCurr);

%Limit vel to 20% search range:
vMax = 0.2;

sV = sign(vel);
rng = repmat((prmBnds(:,2) - prmBnds(:,1))', [szCf(1), 1]);
vel = sV.*min(abs(vel), vMax*rng);

%Create new parameters:
cfN = cfCurr + vel;

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
