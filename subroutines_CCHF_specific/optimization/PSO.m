% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function coefN = pso(coef, fitScore, prmBnds, iter)

%SEE DISCUSSION AT END OF SCRIPT REGARDING PSO OPTIMIZATION.

persistent vel
if iter == 1
   vel = zeros(size(squeeze(coef(1,:,:))));
end

%Set learning rate:
lrnRt = 0.5; %According to source pasted below, common value for scalar is 2
   

%Determine number of Monte Carlo runs (based on number of parameters)
nParam = numel(prmBnds(:,1));

nRCarlo = n_monte_carlo(nParam);

nPop = numel(fitScore(1,:));
nItCarlo = round(nRCarlo/nPop);

if iter == 1
   disp([char(10) 'Based on the number of parameters being optimized (' ...
       num2str(nParam) '), optimization will proceed with ' ...
       num2str(nItCarlo) ' iterations of Monte Carlo simulation prior to PSO.' char(10)]); 
end


%Generate next group of coefficients:
if iter <= nItCarlo %For first couple generations, use Monte Carlo
    coefN = Monte_Carlo(prmBnds, numel(fitScore(1,:)));
else %Regular PSO optimization:
    [coefN, vel] = pso_update(coef(1:iter,:,:), fitScore(1:iter,:), prmBnds, lrnRt, vel);
   
    %Scale velocity "inertia":
    %Method from Y.ShiandR.C.Eberhart,?Amodifiedparticleswarmoptimizer,?inProc. IEEE World Congr. Comput. Intell., 1998, pp. 69?73
    omMin = 0.4;
 	omMax = 0.9;
    om = omMax - (omMax - omMin)*min(iter/numel(coef(:,1,1)), 1);
    vel = om*vel;
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
