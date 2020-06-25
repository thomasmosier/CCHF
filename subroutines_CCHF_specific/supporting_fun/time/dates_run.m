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

function sMeta = dates_run(sMeta, cal, varargin)



% dateRef = [0,1,1,1,1]; %year, month, day, hour, minute


blCellConvert = 0;
if ~iscell(sMeta.dateStart) && ~iscell(sMeta.dateStart)
    dateLpStart = {sMeta.dateStart};
    dateLpEnd   = {  sMeta.dateEnd};
    
    blCellConvert = 1;
elseif iscell(sMeta.dateStart) && iscell(sMeta.dateStart)
    dateLpStart = sMeta.dateStart;
    dateLpEnd   =   sMeta.dateEnd;
else
    error('datesRun:DiffDateFormat','The input date fields have different formats. This has not been programmed for.')
end

dateRun = cell(numel(dateLpStart), 1);

%Loop over each set of input dates
for ii = 1 : numel(dateRun(:))
    %Adjust start and end date vector to have same temporal resolution as timestep
    if regexpbl(sMeta.dt,{'daily','day'})
        if length(dateLpStart{ii}) == 2
            dateLpStart{ii} = [dateLpStart{ii}, 1];
        elseif length(dateLpStart{ii}) > 3
            dateLpStart{ii} = dateLpStart{ii}(1:3);
        elseif length(dateLpStart{ii}) == 3

        else
            error('dates_run:shortDate',['The start date has a length that '...
                'has not been programmed for.']);
        end

        if length(dateLpEnd{ii}) == 2
            dateLpEnd{ii} = [dateLpEnd{ii}, eomday(dateLpEnd{ii}(1), dateLpEnd{ii}(2))];
        elseif length(dateLpEnd{ii}) > 3
            dateLpEnd{ii} = dateLpEnd{ii}(1:3);
        elseif length(dateLpStart{ii}) == 3

        else
            error('dates_run:shortDate',['The end date has a length that '...
                'has not been programmed for.']);
        end
    elseif regexpbl(sMeta.dt,'month')
        if length(dateLpStart{ii}) == 1
            dateLpStart{ii} = [dateLpStart{ii}, 7];
        elseif length(dateLpStart{ii}) == 2

        elseif length(dateLpStart{ii}) > 2
            dateLpStart{ii} = dateLpStart{ii}(1:2);
        else
            error('dates_run:shortDate',['The start date has a length that'...
                ' has not been programmed for.']);
        end

        if length(dateLpEnd{ii}) == 1
            dateLpEnd{ii} = [dateLpEnd{ii}, 6];
        elseif length(dateLpStart{ii}) == 2

        elseif length(dateLpEnd{ii}) > 2
            dateLpEnd{ii} = dateLpEnd{ii}(1:2);
        else
            error('dates_run:shortDate',['The end date has a length that'...
                ' has not been programmed for.']);
        end
    elseif regexpbl(sMeta.dt,'hour')
        if length(dateLpStart{ii}) == 3
            dateLpStart{ii} = [dateLpStart{ii}, 1];
        elseif length(dateLpStart{ii}) == 2
            dateLpStart{ii} = [dateLpStart{ii}, 1, 1];
        elseif length(dateLpStart{ii}) == 4

        elseif length(dateLpStart{ii}) > 4
            dateLpStart{ii} = dateLpStart{ii}(1:4);
        else
            error('dates_run:shortDate',['The start date has a length that '...
                'has not been programmed for.']);
        end

        if length(dateLpEnd{ii}) == 3
            dateLpEnd{ii} = [dateLpEnd{ii}, 24];
        elseif length(dateLpEnd{ii}) == 2
            dateLpEnd{ii} = [dateLpEnd{ii}, eomday(dateLpEnd{ii}(1), dateLpEnd{ii}(2)), 24];
        elseif length(dateLpStart{ii}) == 4

        elseif length(dateLpEnd{ii}) > 4
            dateLpEnd{ii} = dateLpEnd{ii}(1:3);
        else
            error('dates_run:shortDate',['The end date has a length that '...
                'has not been programmed for.']);
        end
    end

    %Add spin to time:
    if ~isempty(varargin(:)) && regexpbl(varargin{1},'spin')
        indSpin = regexpi(sMeta.spinup,'[0-9]');

        if isempty(indSpin)
            error('dates_run:spinupNumber',['No numbers detected in spinup string' char(39) sMeta.spinup char(39)]);
        elseif any(diff(indSpin)> 1)
            error('dates_run:spinupNumber',['Multiple numbers detected in spinup string' char(39) sMeta.spinup char(39)]);
        else
            nSpin = str2double(sMeta.spinup(indSpin(1):indSpin(end)));
        end

        if  regexpbl(sMeta.spinup,'month')
            if length(dateLpStart{ii}) == 2
                dateSpinStart = dateLpStart{ii} - [0, nSpin];
            elseif length(dateLpStart{ii}) == 3
                dateSpinStart = dateLpStart{ii} - [0, nSpin, 0];
            elseif length(dateLpStart{ii}) == 4
                dateSpinStart = dateLpStart{ii} - [0, nSpin, 0, 0];
            else
                error('dates_run:spinupMonth',[length(dateLpStart{ii}) ' for time unit of month is an unrecognized length.']);
            end

            %Correct for negative months:
            if dateSpinStart(2) < 1
               yrExtraFrac = -dateSpinStart(2)/12; 
               dateSpinStart(1) = dateSpinStart(1) - ceil(yrExtraFrac);
               dateSpinStart(2) = mod(dateSpinStart(2)-1,12)+1;
            end
        elseif regexpbl(sMeta.spinup,'year')
            if length(dateLpStart{ii}) == 1
                dateSpinStart = dateLpStart{ii} - nSpin;
            elseif length(dateLpStart{ii}) == 2
                dateSpinStart = dateLpStart{ii} - [nSpin, 0];
            elseif length(dateLpStart{ii}) == 3
                dateSpinStart = dateLpStart{ii} - [nSpin, 0, 0];
            elseif length(dateLpStart{ii}) == 4
                dateSpinStart = dateLpStart{ii} - [nSpin, 0, 0, 0];
            else
                error('dates_run:spinupMonth',[length(dateLpStart{ii}) ' for time unit of year is an unrecognized length.']);
            end
        elseif regexpbl(sMeta.spinup,{'day','daily'})
            if length(dateLpStart{ii}) == 3
                dateSpinStart = dateLpStart{ii} - [sMeta.spinup, 0];
            else
                error('dates_run:spinupMonth',[length(dateLpStart{ii}) ' for time unit of year is an unrecognized length.']);
            end

            if dateSpinStart(3) < 1 
                if eomday(dateLpStart{ii}(1), mod(dateLpStart{ii}(2)-1,12)+1) - dateSpinStart(3) > 0
                    dateLpStart{ii}(2) = dateLpStart{ii}(2) - 1;
                    dateLpStart{ii}(3) = eomday(dateLpStart{ii}(1), mod(dateLpStart{ii}(2)-1,12)+1) - dateSpinStart(3);
                else
                    error('dates_run:spinupDay','The spinup period is in units of days but greater than one month.  This has not been programmed for.')
                end
            end
        else
           error('dates_run:spinupUnit',[sMeta.spinup ' is an unrecognized time unit.']);
        end

    else
        dateSpinStart = dateLpStart{ii};
    end
    
    dateRun{ii} = date_vec_fill(dateSpinStart,dateLpEnd{ii},cal);
end


%Convert output to input format:
if blCellConvert == 1
    sMeta.dateRun = dateRun{1};
else
    sMeta.dateRun = dateRun;
end
