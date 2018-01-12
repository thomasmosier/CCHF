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

function dataOut = interp_month2day(datesMonth, dataIn, datesDays, method, blNorm)


daysInMonthIn = nan(numel(datesMonth(:,1)),1);
for ii = 1 : numel(daysInMonthIn)
   daysInMonthIn(ii) = eomday(datesMonth(ii,1), datesMonth(ii,2));
end
daysIn = [0; cumsum(daysInMonthIn(1:end-1))] + 0.5*daysInMonthIn;

daysInMonthOut = nan(numel(datesDays(:,1)),1);
cumDays = 0;
for ii = 1 : numel(daysInMonthOut)
    if ii ~= 1 && diff([datesDays(ii,2),datesDays(ii-1,2)]) ~= 0
        cumDays = cumDays + eomday(datesDays(ii-1,1), datesDays(ii-1,2));
    end
   daysInMonthOut(ii) = cumDays + datesDays(ii,3);
end

compareFirst = datesDays(1,:) < [datesMonth(1,:), 0.5*eomday(datesMonth(1,1), datesMonth(1,2))];

if isequal(compareFirst,[0,0,1])
    daysOffset = 0;
else
    error('interp_month2day:diffStartDates', ['The input data date and '...
        'interp2 date are not the same.  This has not been programmed for.']);
end

daysIn = daysIn + daysOffset;

if ndims(dataIn) == 2 
    if any(size(dataIn) == 1)
        if sum(~isnan(dataIn)) > 2
            dataOut = interp1(daysIn, dataIn, daysInMonthOut, method);
        else
            dataOut = nan(numel(daysInMonthOut),1);
        end

        %Normalize if specified:
        if blNorm == 1
            for ii = 1 : numel(dataIn)
                indCurr = ismember(datesDays(:,1:2), datesMonth(ii,:),'rows');

                scale = dataIn(ii)/nansum(dataOut(indCurr));
                scale(isnan(scale)) = 0;
                dataOut(indCurr) = scale*dataOut(indCurr);
            end
        elseif blNorm ~= 0
            error('interp_month2day:normalize',['Normalization value is ' ...
                num2str(blNorm) ' but must be either 0 or 1.']);
        end
    else
        error('interp_month2day:inputMatrix','The input data is a matrix with no singleton dimensions, which has not been programmed for.')
    end
elseif ndims(dataIn) == 3
    nLat = numel(dataIn(1,:,1));
    nLon = numel(dataIn(1,1,:));
    dataOut = nan([numel(daysInMonthOut), nLat, nLon], 'single');
    
    for ii = 1 : nLat
        for jj = 1 : nLon
            if sum(~isnan(dataIn)) > 2
                dataOut(:,ii,jj) = interp1(daysIn, squeeze(dataIn(:,ii,jj)), daysInMonthOut, method);
            else
                dataOut(:,ii,jj) = nan(numel(daysInMonthOut),1);
            end
        end
    end

    %Normalize if specified:
    if blNorm == 1
        for ii = 1 : numel(dataIn(:,1,1))
            indCurr = find(ismember(datesDays(:,1:2), datesMonth(ii,:),'rows') == 1);

            scale = dataIn(ii,:,:)./nansum(dataOut(indCurr,:,:), 1);
            scale(isnan(scale)) = 0;
            for jj = 1 : numel(indCurr)
                dataOut(indCurr(jj),:,:) = squeeze(scale).*squeeze(dataOut(indCurr(jj),:,:));
            end
        end
    elseif blNorm ~= 0
        error('interp_month2day:normalize',['Normalization value is ' ...
            num2str(blNorm) ' but must be either 0 or 1.']);
    end
else
   error('interp_month2day:unknownDim', ['The current data has ' ...
       num2str(ndims(dataIn)) ' dimensions, which has not been programmed for.']); 
end
    

% plot(daysInMonthOut, dataOut, '.', daysIn, dataIn, '*');