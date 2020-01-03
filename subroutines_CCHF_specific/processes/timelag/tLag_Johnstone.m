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

% This may not work well for large spatial domains

sztLag = size(sHydro.dl);


%Determine if array exists or needs to be calculated:
if ~isfield(sLand, 'tlag')
    %Create path for saving/loading tlag (this matters because for large
    %domains can take many hours to calculate)
    dirTlag = fullfile(sMeta.foldstorage, 'tLag_Johnstone');
    if ~exist(dirTlag, 'dir')
       mkdir(dirTlag); 
    end
    pathTlag = fullfile(dirTlag, [sMeta.region{sMeta.siteCurr} '_timelag_Johnstone_pwr=' num2str(round2(a,3)), '.mat']);

    varTLag = 'tlag';

    %Boolean value to indicate that time-lag array must be calculated
    blCalc = 0;
    
    if ~isempty(pathTlag) && exist(pathTlag, 'file')
        SLd = load(pathTlag);
        
        if isfield(SLd, varTLag)
            sLand.tlag = SLd.(varTLag);
        else
            blCalc = 1;
        end
    else
        blCalc = 1;
    end
    
    %Calculate array if needed.
    if blCalc == 1
        sLand.tlag = sparse(sztLag(1),sztLag(1));
        iTLag = find(sHydro.dl ~= 0);
        tic
        sLand.tlag(iTLag) = 10^a*3600*sqrt(sHydro.dl(iTLag) ./ sqrt(abs(sHydro.slope(iTLag))));

        deltaT = toc;
        if deltaT > 20*60
            disp(['It took ' num2str(toc/60) ' minutes to define time-lag array.']);
        end
        
        tlag = sLand.tlag;
        if ~isempty(pathTlag)
            save(pathTlag, 'tlag')
        end
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
end
    