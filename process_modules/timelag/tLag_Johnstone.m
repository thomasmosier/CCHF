function varargout = tLag_Johnstone(sHydro, varargin)
%Follows lag formulation from:
%Bedient, P. B., & Huber, W. C. (1988). Hydrology and floodplain analysis.
%pg. 129.
%Which cites:
%Johnstone, D., & Cross, W. P. (1949). Elements of applied hydrology.


global sLand


if isempty(varargin(:))
    varargout{1} = cell(0,6);
    varargout{1} = cat(1, varargout{1}, {'timelag_pwr', -2.5, 1, -0.5, 'tLag_Johnstone','routing'}); %Unitless scalar
%     argout = cell(1,6);
%     argout(1,:) = {'timelag_pwr', 0.5, 1.5, 0.9561, 'tLag_Clark', 'routing'}; %Unitless scalar
    return
else
    sMeta = varargin{1};
    a = find_att(sMeta.coef,'timelag_pwr'); %Unitless scalar
end

%MAKE tLag same size as DEM (can't handle large sparse matrice because of
%memory issues:
% sLand.tlagFdr = full(reshape(sum(sHydro.fdr.*sLand.tlag, 2), size(sHydro.dem)));

% This doesn't work for large spatial domains:
sztLag = size(sHydro.dl);


% nm = [sMeta.region '-Johnstone_tLag_array-' num2str(round2(a,3))];
% pathTLag = fullfile(sMeta.rtDir,[nm '.asc']);



sLand.tlag = sparse(sztLag(1),sztLag(1));
iTLag = find(sHydro.dl ~= 0);
tic
sLand.tlag(iTLag) = 10^a*3600*sqrt(sHydro.dl(iTLag) ./ sqrt(abs(sHydro.slope(iTLag))));

deltaT = toc;
if deltaT > 20*60
    disp(['It took ' num2str(toc/60) ' minutes to define time-lag array.']);
end


% if exist(pathTLag,'file') %Save time by reading time lag array from file (matters in case of large spatial domains
%     [dataTemp,~,~] = read_ESRI(pathTLag);
%     if isequal(size(dataTemp),sztLag)
%         sLand.tlag = sparse(dataTemp);
%     else
%         error('tLag_Johnstone:sizeLoaded','The size of the loaded time lag array is incorrect.');
%     end
%     clear('dataTemp');
% else
%     sLand.tlag = sparse(sztLag(1),sztLag(1));
%     iTLag = find(sHydro.dl ~= 0);
%     tic
%     sLand.tlag(iTLag) = 10^a*3600*sqrt(sHydro.dl(iTLag) ./ sqrt(abs(sHydro.slope(iTLag))));
%     
% %     %Write timelag to save time during future runs.
% %     write_ESRI_v4(sLand.tlag, nan(6,1), pathTLag, 0);
%     deltaT = toc;
%     if deltaT > 20*60
%         disp(['It took ' num2str(toc/60) ' minutes to define time-lag array.']);
%     end
% end


% ntLag = numel(sHydro.dl);
% if ntLag > 10^6
%    step = floor(10^7/sztLag(1));
%    nStep = ceil(sztLag(2)/step);
%    cntr = 1;
%    
% 	for ii = 1 : nStep
%         if ii == nStep
%             sLand.tlag(:,cntr:end) ...
%                 = 10^a*3600*sqrt(sHydro.dl(:,cntr:end) ...
%                 ./ sqrt(abs(sHydro.slope(:,cntr:end))));
%         else
%             sLand.tlag(:,cntr:cntr+step-1) ...
%                 = 10^a*3600*sqrt(sHydro.dl(:,cntr:cntr+step-1) ...
%                 ./ sqrt(abs(sHydro.slope(:,cntr:cntr+step-1))));
%         end
%         
%         cntr = cntr + step;
% 	end
% else
% 	sLand.tlag(iTLag) = 10^a*3600*sqrt(sHydro.dl(iTLag) ./ sqrt(abs(sHydro.slope(iTLag))));  %(units = seconds)
% end




%Trim t-lag entries:
%Set Nan's to 0:
sLand.tlag(isnan(sLand.tlag)) = 0;
sLand.tlag(sLand.tlag == Inf) = max2d(sLand.tlag);

    