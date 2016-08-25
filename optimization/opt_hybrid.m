function coefN = opt_hybrid(coef, fitScore, prmBnds, iter)

%SEE DISCUSSION AT END OF SCRIPT REGARDING PSO OPTIMIZATION.

%Determine number of Monte Carlo runs (based on number of parameters)
nParam = numel(prmBnds(:,1));
if nParam <= 4
    nRCarlo = 150;
elseif nParam <= 7
    nRCarlo = 300;
elseif nParam > 7 && nParam < 10
    nRCarlo = 500;
else
    nRCarlo = 800;
end

nPop = numel(fitScore(1,:));
nItCarlo = round(nRCarlo/nPop);

if iter == 1
   disp(['Based on the number of parameters being optimized (' ...
       num2str(nParam) '), optimization will proceed with ' ...
       num2str(nItCarlo) ' iterations of Monte Carlo simulation prior to PSO.']); 
end

%Scalar that determines speed of movement in PSO 
scalar = 0.4; %According to source pasted below, common value for scalar is 2
    

%Keep only current and past generation from inputs:
coef = coef(1:iter,:,:);
fitScore = fitScore(1:iter,:);

valMostFit = round2(min(fitScore,[],2),2);
%Determine how many iterations the fitness score has been at the current
%value:
[~, indGenFit] = min(valMostFit);


%Generate next group of coefficients:
if iter <= nItCarlo %For first couple generations, use Monte Carlo
    disp('Monte Carlo Sampling used to explore global parameter space.')
    coefN = Monte_Carlo(prmBnds,numel(fitScore(1,:)));
elseif iter > 5 + nItCarlo && indGenFit + 7 < iter && indGenFit > nItCarlo
    disp('This iteration was optimized with linear-sensitivity to explore local parameter space.')
    fDiv = 0.4*rand;
    
    
    coefN = Monte_Carlo(prmBnds,numel(fitScore(1,:)));
    nPop = numel(fitScore(1,:));
    [~, pLBest] = min(fitScore(iter,:));
    if nPop >= 2*nParam
        cntr = 0;
        for jj = 1 : 2
            for ii = 1 : nParam
                cntr = cntr + 1;
                
                coefN(cntr,:) = squeeze(coef(iter, pLBest,:));
                if jj == 1
                    coefN(cntr,ii) = coefN(cntr,ii) + fDiv*(prmBnds(ii,2)-prmBnds(ii,1));
                else
                    coefN(cntr,ii) = coefN(cntr,ii) - fDiv*(prmBnds(ii,2)-prmBnds(ii,1));
                end
            end
        end
    elseif nPop >= nParam
        cntr = 0;
        for ii = 1 : nParam
            cntr = cntr + 1;

            coefN(cntr,:) = squeeze(coef(iter, pLBest,:));
            if rand > 0.5
                coefN(cntr,ii) = coefN(cntr,ii) + fDiv*(prmBnds(ii,2)-prmBnds(ii,1));
            else
                coefN(cntr,ii) = coefN(cntr,ii) - fDiv*(prmBnds(ii,2)-prmBnds(ii,1));
            end
        end
    else
        [~, ord] = sort(rand(nParam,1));
        cntr = 0;
        for ii = 1 : nPop
            cntr = cntr + 1;

            coefN(cntr,:) = squeeze(coef(iter, pLBest,:));
            if rand > 0.5
                coefN(cntr,ord(ii)) = coefN(cntr,ord(ii)) + fDiv*(prmBnds(ord(ii),2)-prmBnds(ord(ii),1));
            else
                coefN(cntr,ord(ii)) = coefN(cntr,ord(ii)) - fDiv*(prmBnds(ord(ii),2)-prmBnds(ord(ii),1));
            end
        end
    end
    
    %Fittest parameter set from current generation evolves in PSO manner:
    [~, pLBest] = min(fitScore(iter,:));
    [~, iGBest, pGBest] = min2d(fitScore);
    
    pv = scalar*(squeeze(coef(iGBest,pGBest,:)) - squeeze(coef(iter,pLBest,:)));
    
%     pv = scalar*(squeeze(coef(rBest,cBest,:)) - squeeze(coef(iter,iGBest,:))) ...
%         ./(prmBnds(:,2)-prmBnds(:,1));
%     pv(pv < -1) = -1;
%     pv(pv > 1) = 1;
    
    coefN(pLBest,:) = squeeze(coef(iter,pLBest,:)) + pv;
elseif std(fitScore(iter,:)) < 0.3 %If fitness scores converge too much, add random diversity
    disp('This iteration was optimized with Monte Carlo sampling to add diversity.')
    coefN = Monte_Carlo(prmBnds,numel(fitScore(1,:)));
    
    %Fittest parameter set from current generation evolves in PSO manner:
    [~, pLBest] = min(fitScore(iter,:));
    [~, iGBest, pGBest] = min2d(fitScore);
    
    pv = scalar*(squeeze(coef(iGBest,pGBest,:)) - squeeze(coef(iter,pLBest,:)));
    
    coefN(pLBest,:) = squeeze(coef(iter,pLBest,:)) + pv;
else %Regular PSO optimization:
    disp('This iteration was optimized with PSO algorithm.')
    %Identify best parameter set:
    [~, iPBest] = min(fitScore,[],1);
    [~, pLBest] = min(fitScore(iter,:));

    sz = [numel(coef(1,:,1)), numel(coef(1,1,:))];
    cfPBest = nan(sz);
    for ii = 1 : sz(1)
        cfPBest(ii,:) = coef(iPBest(ii),ii,:);
    end
    
    %Calculate velocity for each parameter sets:
    pv = scalar*((cfPBest - squeeze(coef(iter,:,:))) ...
        + (ones(numel(fitScore(1,:)),1)*squeeze(coef(iter,pLBest,:))' - squeeze(coef(iter,:,:))));
%     v = scalar*((cfPBest - squeeze(coef(iter,:,:))) ...
%         + (ones(numel(fitScore(1,:)),1)*squeeze(coef(iter,iGBest,:))' - squeeze(coef(iter,:,:))));
%     pv = v./(ones(numel(fitScore(1,:)),1)*(prmBnds(:,2)-prmBnds(:,1))');
%     pv(pv < -1) = -1;
%     pv(pv > 1) = 1;

    %Generate new parameter sets:
    coefN = squeeze(coef(iter,:,:)) + pv;
    for ii = 1 : sz(2)
       indNeg = find(coefN(:,ii) < prmBnds(ii,1));
       indPos = find(coefN(:,ii) > prmBnds(ii,2));
       if ~isempty(indNeg)
           coefN(indNeg,ii) = prmBnds(ii,1);
       end
       if ~isempty(indPos)
           coefN(indPos,ii) = prmBnds(ii,2);
       end
    end
end



%%FROM http://www.swarmintelligence.org/tutorials.php
% 3. The algorithm
% As stated before, PSO simulates the behaviors of bird flocking. Suppose 
%the following scenario: a group of birds are randomly searching food in an 
%area. There is only one piece of food in the area being searched. All the 
%birds do not know where the food is. But they know how far the food is in 
%each iteration. So what's the best strategy to find the food? The 
%effective one is to follow the bird which is nearest to the food. 

% PSO learned from the scenario and used it to solve the optimization 
%problems. In PSO, each single solution is a "bird" in the search space. We 
%call it "particle". All of particles have fitness values which are 
%evaluated by the fitness function to be optimized, and have velocities 
%which direct the flying of the particles. The particles fly through the 
%problem space by following the current optimum particles. 
% 
% PSO is initialized with a group of random particles (solutions) and then 
%searches for optima by updating generations. In every iteration, each 
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
% While maximum iterations or minimum error criteria is not attained
% 
% Particles' velocities on each dimension are clamped to a maximum velocity 
%Vmax. If the sum of accelerations would cause the velocity on that 
%dimension to exceed Vmax, which is a parameter specified by the user. Then 
%the velocity on that dimension is limited to Vmax.
