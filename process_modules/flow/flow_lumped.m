function flow_lumped(sHydro, varargin)
%sLand.flow has units of m^3 / time unit used if sLand.mrro has units of m.


global sLand

%If not variable input argument, return parameters and do not proceed.
%Clark time lag:
if isempty(varargin(:))
    return
else
    sMeta = varargin{1};
end

if sum(diag(sHydro.fdr)) > 0
   warning('flow_Clark:endoheroicFlow',['Some cells flow to themselves.  This may '...
       'cause an infinite loop or may not be handled properly.']); 
end


%Find dt (seconds):
dt = time2sec(1,sMeta.dt,sMeta.dateCurr); %Used for distributing previous time-steps streamflow (units = seconds)


if ~isfield(sLand,'totFdrR') || ~isfield(sLand,'totFdrF') || ~regexpbl(sMeta.module,'time-Johnstone')
    szRun = size(sLand.mrro);
    cumTime  = zeros(szRun ,'single'); %(units = seconds)
    fdrSteps = zeros(szRun,'single');
    
    [sLand.totFdrR, fdrSteps, cumTime] = fdr_cum(sHydro, sLand.tlag, fdrSteps, cumTime,   dt);
    [sLand.totFdrF,        ~,       ~] = fdr_cum(sHydro, sLand.tlag, fdrSteps, cumTime, 2*dt);
end


% if isequal(sMeta.dateCurr,[2001, 7,30])
%     keyboard
% end
if sMeta.indCurr == 1
    %Add runoff to flow for current time step:
    sLand.flow(:) = (speye(size(sLand.totFdrR)) + sLand.totFdrR)'*double(reshape(sLand.mrro,[],1))/dt;
else
    sLand.flow(:) = (speye(size(sLand.totFdrR)) + sLand.totFdrR)'*double(reshape(sLand.mrro,[],1))/dt ...
        + sLand.totFdrF'*double(sLand.flow(:));
end
%In FDR: row is cell of origin and col is downstream cell



%%ROUTE FLOW FROM PREVIOUS TIME-STEP
%
%ONLY ROUTE FLOW FROM CELLS HIGHER UP THAN WOULD HAVE BEEN INCLUDED IN
%ROUTING CURRENT TIMESTEP (uses fdrStep to accomplish this and does not
%reset 'ii')
% if sLand.indCurr ~= 1
% %     %Create indexes for flow from previous time-step:
% %     indFlowFrom = (sLand.indCurr*ones(nGrid,1) - 1) + (rFlow-1)*szFlow(1) + (cFlow-1)*szFlow(2)*szFlow(1);
% %     
%     %Add previous timestep flow to current timestep:
%     sLand.flow(:) = sLand.flow(:) + sLand.totFdrF'*sLand.flow(:);
% end

%Flow from multiple previous timesteps:
% if sLand.indCurr ~= 1
%     totFdr = zeros(size(sHydro.fdr));
% 
%     %First index is that of cell upstream of runoff routing from past
%     %time-step:
%     ii = min2d(fdrSteps) + 1;
%     
%     flagFl = 1; tt = 0; %Don't include flow at each cell from previous time step (but route upstream flow)
%     while flagFl == 1 && ii < nGrid %Loop over time step that the flow is flowing from
%         tt = tt + 1;
%         %Create indexes for flow from previous time-step:
%         indFlowFrom = (sLand.indCurr*ones(nGrid,1) - tt) + (rFlow-1)*szFlow(1) + (cFlow-1)*szFlow(2)*szFlow(1);
%         
%         while flagFl == 1 && ii < nGrid %Create a fdr grid for the current set of flow-from -> flow-to points
%             currFdr = full(sHydro.fdr^ii);   
% 
%             indFdrSkip = find(fdrSteps + 1 > ii);
%             indCurrFdr = setdiff(indRun, indFdrSkip);
%             %Add time lag to cumulative time array
%             tTravCurr = diag(currFdr*sLand.tlag,0);
% 
%             cumTime(indCurrFdr) = cumTime(indCurrFdr) + tTravCurr(indCurrFdr); %(units = seconds)
% 
%     %         indRem = find(cumTime <= dt | cumTime > 2*dt); %Ensure 
% 
%             currFdr(cumTime > tt*dt,:) = 0;
%             currFdr(indFdrSkip,:) = 0;
%             ii = ii + 1;
% 
%             totFdr = totFdr + currFdr;
%             fdrSteps(indCurrFdr) = fdrSteps(indCurrFdr) + 1;
% 
%             %Exit when no more water to route or tiem has been exceeded at all
%                 %points:
%             if all(all(currFdr == 0)) || all(all(cumTime > tt*dt))
%                flagFl = 0; 
%             end
%         end
% 
%         %Add previous timestep flow to current timestep:
%         sLand.flow(indFlowTo) = sLand.flow(indFlowTo) + totFdr'*sLand.flow(indFlowFrom);
%     end
% end