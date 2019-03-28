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

function varargout = heat_SETI_Mosier(varargin)
%See reference:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, flux_long_simple());
%     varargout{1} = cat(1,varargout{1}, flux_long_Mosier(sHydro));
%     varargout{1} = cat(1,varargout{1}, flux_long_Deacon(sHydro));
%     varargout{1} = cat(1,varargout{1}, flux_precip());
	varargout{1} = cat(1,varargout{1}, flux_conduct_snow_linear()); 
%     varargout{1} = cat(1,varargout{1}, flux_conduct_c()); 
    varargout{1} = cat(1,varargout{1}, {'Watt_per_deg',  0, 35, 10, 'heat_SETI_Mosier','cryo'});
    varargout{1} = cat(1,varargout{1}, ['tsn', cell(1,5)]);
        
    return
else
    sMeta = varargin{1};
    wattperdeg = find_att(varargin{1}.coef,'Watt_per_deg'); 
%     lrf = find_att(varargin{1}.coef,'lw_pwr'); 
%     srf = find_att(varargin{1}.coef,'SW_pwr'); 
end

% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
%Units of J-m^{-2}-hr-mm^{-1}-s^{-1} := W-m^{-2}-hr-mm^{-1} (because [dis] = mm/hr/deg 

minTemp = find_att(sMeta.global,'snow_temp_min');

%Initialize skin temperature for snow and ice
if ~isfield(sCryo,'tssn') 
   sCryo.tssn = zeros(size(sCryo.snw),'single'); 
end
if ~isfield(sCryo,'tsic') 
   sCryo.tsic = zeros(size(sCryo.snw),'single'); 
end
if ~isfield(sCryo,'tsn') 
   sCryo.tsn = zeros(size(sCryo.snw),'single'); 
end
if ~isfield(sCryo,'snlq') 
   sCryo.snlq = zeros(size(sCryo.snw),'single'); 
end


%Find top of the atmosphere radiation:
indRSDT = find(ismember(sAtm.datersdt(:,2:end), sMeta.dateCurr(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end
    
%Shortwave melt energy (independent of surface temperature):
sCryo.hfrs  = (1-sCryo.snalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));
sCryo.hfrsi = (1-sCryo.icalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));


%%IMPLEMENT BISECTION METHOD TO ITERATIVELY SOLVE FOR SURFACE TEMPERATURE
%positive indicates melt occuring

% if isequal(sMeta.dateCurr(1:2),[2000, 1])
%     keyboard
% end


%*Could implement Brent Method (https://en.wikipedia.org/wiki/Brent%27s_method)
% tic
sCryo.tssn = roots_bisect(minTemp,'s');
sCryo.tsic = roots_bisect(minTemp,'i');
% toc

% tic
% sCryo.tsis = roots_Newton_Raphson(minTemp);
% toc

%Set minimum snow surface temperature:
sCryo.tssn( sCryo.tssn < minTemp) = minTemp; 
sCryo.tsic( sCryo.tsic < minTemp) = minTemp; 
%Set snow temperature to nan where there is no snow:
sCryo.tssn(sCryo.snw == 0) = nan;
sCryo.tsic(sCryo.icx == 0) = nan;
%Set maximum snow temperature to 0
sCryo.tssn( sCryo.tssn > 0) = 0;
sCryo.tsic( sCryo.tsic > 0) = 0;

%Set surface temp to middle value:
sCryo.hfnet = SETI_heatflux(sCryo.tssn,'s','norm','global');
sCryo.hfneti = SETI_heatflux(sCryo.tsic,'i','norm','global');



%DEFINE NET HEAT FUNCTION
    function ht = SETI_heatflux(tsis, surLyr, varargin)
    %     flux_precip(varargin{1});

        if isempty(varargin(:)) || ~regexpbl(varargin{1}, 'deriv') %Normal evaluation
            hft = wattperdeg*(squeeze(sAtm.tas(sAtm.indtas,:,:)) - tsis);
            hfrl = flux_long_simple(tsis, sMeta);
            
            if regexp(surLyr,'s')
                hfsnc = flux_conduct_snow_linear(tsis, sMeta);
                
                ht = sCryo.hfrs + hfrl + hft + hfsnc;
            elseif regexp(surLyr,'i') %Different shortwave radiation and conduction for glaciers
                ht = sCryo.hfrsi + hfrl + hft;
        %         flux_glacier_conduct(varargin{1});
        %         ht = sCryo.hfrsi + sCryo.hfrl + sCryo.hft + sCryo.hficc;
            else
                error('SETI_Mosier:unknownSurfType', [surLyr ' is an unknown input.']);
            end  
        elseif ~isempty(varargin(:)) && regexpbl(varargin{1}, 'deriv') %Find derivative of heatflux
            dhfrl = flux_long_simple(tsis, sMeta, 'deriv');
            dhft = -wattperdeg;
        
            if regexp(surLyr,'s')
                dhfsnc = flux_conduct_snow_linear(tsis, sMeta, 'deriv');
                ht = dhfrl + dhft + dhfsnc;
            elseif regexp(surLyr,'i') %Different shortwave radiation and conduction for glaciers
                ht = dhfrl + dhft;
        %         flux_glacier_conduct(varargin{1});
        %         ht = sCryo.hfrsi + sCryo.hfrl + sCryo.hft + sCryo.hficc;
            else
                error('SETI_Mosier:unknownSurfType', [surLyr ' is an unknown input.']);
            end
        else
            error('SETI_Mosier:unknownMode', [varargin{1} ' is an unknown mode.']);
        end
        
        %If option selected, save fields to global variable:
        if numel(varargin(:)) > 1 && regexpbl(varargin{2}, 'glob')
            sCryo.hfrl = hfrl;
            sCryo.hft = hft;
            
            if regexp(surLyr,'s')
                sCryo.hfsnc = hfsnc;
            end
        end
            
    end

%FUNCTIONS TO ITERATIVELY SOLVE FOR SURFACE TEMPERATURE:

%Bisection method (no derivates needed):
    function tsnNew = roots_bisect(minTemp, surLyr)
        
        nLp = 100;
        tol = 0.1;
        cntr = 0;
        indAll = (1:numel(sCryo.snw));

        %Define endpoints:
        %Temp High
        % tsnH = 40*ones(size(sCryo.snw),'single');
        tsnH = 5*ones(size(sCryo.snw),'single');
        %Temp Low
        tsnL = minTemp*ones(size(sCryo.snw),'single');


        %Iteratively solve for surface temperature USING BISECTION METHOD:
        while cntr <= nLp
            cntr = cntr + 1;

            tsnNew = 0.5*(tsnL+tsnH);

            heatM = SETI_heatflux(tsnNew, surLyr);
            if all2d(abs(tsnH-tsnNew) < tol | isnan(heatM))
               break 
            end

            indSmSgn = find(heatM.*SETI_heatflux(tsnL, surLyr) > 0);
            indOpSgn = setdiff(indAll,indSmSgn);

            tsnL(indSmSgn) = tsnNew(indSmSgn); 
            tsnH(indOpSgn) = tsnNew(indOpSgn); 
        end
    end

%Find Roots using Newton-Raphson Method:
    function tsisNew = roots_Newton_Raphson(minTemp, surLyr)
        nLp = 1000;
        tol = 0.1;
        cntr = 0;

        tsisOld = zeros(size(sCryo.snw));

        %Iteratively solve for surface temperature:
        while cntr <= nLp
            cntr = cntr + 1;

            ht  = SETI_heatflux(tsisOld, surLyr);
            dht = SETI_heatflux(tsisOld, surLyr, 'derivative');
            
            tsisNew = tsisOld - ht./dht;

            if all2d( abs(tsisNew - tsisOld) < tol | isnan(tsisNew))
               break 
            end

        end
        
        tsisNew(tsisNew < minTemp) = minTemp;
    end
end


% %BRENT METHOD FOR ITERATIVELY SOLVING (INCOMPLETE):
% 
% %Ensure |f(L)| < |f(H)|
% indSwap = find(abs(SETI_heat(tsnL,'s')) < abs(SETI_heat(tsnH,'s')));
%     tsnL(indSwap) = tsnH(indSwap);
%     tsnH(indSwap) = minTemp;
%     
% tsnM = tsnL;
% heatM = SETI_heat(tsnM,'s');
% flagM = 1;
% while cntr <= nLp && ~all2d(heatM < tol | isnan(heatM))
%     cntr = cntr + 1;
%     heatL = SETI_heat(tsnL,'s');
%     heatH = SETI_heat(tsnH,'s');
%     
%     if
%         
%     end
%     
%     heatM = SETI_heat(tsnM,'s');
% 
%     
%     indSmSgn = find(SETI_heat(tsnM,'s').*SETI_heat(tsnL,'s') > 0);
%     indOpSgn = setdiff(indAll,indSmSgn);
%     
%     tsnL(indSmSgn) = tsnM(indSmSgn); 
%     tsnH(indOpSgn) = tsnM(indOpSgn); 
% end