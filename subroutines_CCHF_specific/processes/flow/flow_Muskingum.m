function flow_Muskingum(sHydro, varargin)
%sLand.flow has units of m^3 / time unit used if sLand.mrro has units of m.


global sLand

%If not variable input argument, return parameters and do not proceed.
%Clark time lag:
if ~isempty(varargin(:))
    sMeta = varargin{1};
end

if sum(diag(sHydro.fdr)) > 0
   warning('flow_Clark:endoheroicFlow',['Some cells flow to themselves.  This may '...
       'cause an infinite loop or may not be handled properly.']); 
end


%Findsize of flow:
szFlow = size(sLand.flow);

%Find dt (seconds):
dt = time2sec(1,sMeta.dt,sMeta.dateCurr); %Used for distributing previous time-steps streamflow (units = seconds)


%%ROUTE RUNOFF from current time step
%NOTE:
%In fdr, row = grid cell water going from and column = grid cell water
%going to
%To route water from upstream, must have (runoff as column vector)*fdr

%Create routing indices:
% indRun = (1:nGrid)';
% cumTime  = zeros(     szFlow(2:3),'single'); %(units = seconds)
% totFdr   = zeros(size(sHydro.fdr),'single');
% fdrSteps = zeros(     szFlow(2:3),'single');
% 
% 
% flagRf = 1; ii = 1;  
% while flagRf == 1 && ii < nGrid
%     currFdr = full(sHydro.fdr^ii);   
%     
%     %Add time lag to cumulative time array
%     cumTime(indRun) = cumTime(indRun) + diag(currFdr(indRun,:)*sLand.tlag(:,indRun),0); %(units = seconds)
%     fdrSteps(indRun) = ii;
%     
%     indRem = find(cumTime(indRun) >= dt);
%     indRun(indRem) = [];
%     currFdr(indRem,:) = 0;
%     ii = ii + 1;
%     
%     totFdr = totFdr + currFdr;
%     
%     %Exit when no more water to route or time has been exceeded at all
%     %points:
%     if all(all(currFdr == 0)) || isempty(indRun)
%        flagRf = 0; 
%     end
% end

if ~isfield(sLand, 'totFdrF')
    [sLand.totFdrF, ~, ~] = fdr_cum(sHydro, sLand.tlag, zeros(szFlow,'single'), zeros(szFlow,'single'), dt);
end



%Calculate current inflow and previous outflows:
if sLand.indCurr == 1
    sLand.inflowP = zeros(size(sLand.mrro), 'single');
    sLand.flowP = zeros(size(sLand.mrro), 'single');
end

if ~isfield(sLand,'tLagFdr')
    if isfield(sLand, 'tLag')
        sLand.tlagFdr = full(reshape(sum(sHydro.fdr.*sLand.tlag, 2), size(sHydro.dem)));
    else
        error('flow_Muskingum:tLagMissing','The field sLand.tlag is not present but needs to be.');
    end
end

% X = 0.15; %Should be between 0 and 0.3.
% c1 = (dt-2*sLand.tlagFdr*X)./(2*sLand.tlagFdr*(1-X)+dt);
% c2 = (dt+2*sLand.tlagFdr*X)./(2*sLand.tlagFdr*(1-X)+dt);
% c3 = (2*sLand.tlagFdr*(1-X)-dt)./(2*sLand.tlagFdr*(1-X)+dt);
% 
% sLand.outflow = c1.*sLand.inflow + c2.*sLand.inflowP + c3.*sLand.outflowP;

[~, indOrd] = sort(sHydro.flowOrder(:));



X = 0.1;
for ii = 1 : numel(sLand.mrro)
    c1 =     (dt-2*sLand.tlagFdr(indOrd(ii))*X)./(2*sLand.tlagFdr(indOrd(ii))*(1-X)+dt);
    c2 =     (dt+2*sLand.tlagFdr(indOrd(ii))*X)./(2*sLand.tlagFdr(indOrd(ii))*(1-X)+dt);
    c3 = (2*sLand.tlagFdr(indOrd(ii))*(1-X)-dt)./(2*sLand.tlagFdr(indOrd(ii))*(1-X)+dt);
    
    if sum(sLand.totFdrF(:,indOrd(ii))) == 0
        %inflow = runoff at current cell:
       inflow = sLand.mrro(indOrd(ii))/dt;
    else
        %inflow = runoff @ current cell + outflow from upstream cell
        inflow = sLand.mrro(indOrd(ii))/dt ...
            + sum(sLand.flow(sHydro.fdr(:,indOrd(ii)) == 1));
    end
    
    if sLand.indCurr == 1
        sLand.flow(indOrd(ii)) = c1*inflow + c2*sLand.inflowP(indOrd(ii));
    else
        sLand.flow(indOrd(ii)) = c1*inflow + c2*sLand.inflowP(indOrd(ii)) + c3*sLand.flowP(indOrd(ii));
    end

    sLand.inflowP(indOrd(ii)) = inflow;
    %Outflow from current cell must be available as inflow for cells lower
    %down the watershed.
end

%Set current flow as previous outflow for next time step: 
sLand.flowP = sLand.flow;

% X = 0.1;
% for ii = 1 : nGrid
%     c1 =     (dt-2*sLand.tlagFdr(indOrd(ii))*X)./(2*sLand.tlagFdr(indOrd(ii))*(1-X)+dt);
%     c2 =     (dt+2*sLand.tlagFdr(indOrd(ii))*X)./(2*sLand.tlagFdr(indOrd(ii))*(1-X)+dt);
%     c3 = (2*sLand.tlagFdr(indOrd(ii))*(1-X)-dt)./(2*sLand.tlagFdr(indOrd(ii))*(1-X)+dt);
%     
%     if sum(totFdr(:,indOrd(ii))) == 0
%         %inflow = runoff at current cell:
%        inflow = sLand.mrro(indOrd(ii)).*sHydro.area(indOrd(ii))/dt;
%     else
%         %inflow = runoff @ current cell + outflow from upstream cell
%         [rUp, cUp] = ind2sub(szFlow(2:3), find(sHydro.fdr(:,indOrd(ii)) == 1)); 
%         indUp = sLand.indCurr + (rUp-1)*szFlow(1) + (cUp-1)*szFlow(2)*szFlow(1);
%         
%         inflow = sLand.mrro(indOrd(ii)).*sHydro.area(indOrd(ii))/dt ...
%             + sum(sLand.flow(indUp));
%     end
%     
%     [rP, cP] = ind2sub(szFlow(2:3), indOrd(ii)); 
%     indTo = sLand.indCurr + (rP-1)*szFlow(1) + (cP-1)*szFlow(2)*szFlow(1);
%     if sLand.indCurr == 1
%         sLand.flow(indTo) = c1*inflow + c2*sLand.inflowP(indOrd(ii));
%     else
%         indP = sLand.indCurr - 1 + (rP-1)*szFlow(1) + (cP-1)*szFlow(2)*szFlow(1);
%         sLand.flow(indTo) = c1*inflow + c2*sLand.inflowP(indOrd(ii)) + c3*sLand.flow(indP);
%     end
% 
%     sLand.inflowP(indOrd(ii)) = inflow;
%     %Outflow from current cell must be available as inflow for cells lower
%     %down the watershed.
% end