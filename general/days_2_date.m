function [varargout] = days_2_date(daysIn, dateRef, strCal)
% Copyright 2013, 2014 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the GlobalClimateData Downscaling Package.
% 
% The GlobalClimateData Downscaling Package is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The Downscaling Package is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


%Date Vectors have the form: [year, month, day].  
%Can handle multiple dates (where each date is it's own row)

%Change strCal into numerical value defined in 'cal_class'
blCal = cal_class(strCal);

%Make row vector if column:
if length(dateRef(:,1)) > 1 && length(dateRef(1,:)) == 1
    dateRef = dateRef';
end
if length(daysIn(1,:)) > 1 && length(daysIn(:,1)) == 1
    daysIn = daysIn';
end

%Check that 'daysIn' is vector and not matrix
if sum(size(daysIn)>1) >1
   error('days_2_date:daysMatrix',[char(39) 'DaysIn' char(39) ' must be vector, not matrix.']);
end
    
%Initialize date output vector:
%If only one output argument and only one daysIn, then output is [year, month, day] 
if nargout <= 1
    nArg = numel(dateRef);
else
    nArg = nargout;
end
varargout = cell(nArg,1);
[varargout{:}] = deal(nan(numel(daysIn),1));


for ii = 1 : numel(daysIn)
    %Initialize cumulative days:
    daysCurr = 0;
    
    %If either 'daysIn' of the month/year of the reference date are NaN, output NaN:
    if sum([isnan(daysIn(ii,:)), isnan(dateRef(1:2)) ]) > 0
        continue
    else
        for jj = 1 : numel(dateRef)
            varargout{jj}(ii) = dateRef(jj);
        end
    end
    
    %Reset current ouput day to 1 (from daysRef) so that calculation 
    % of days in year/month happens from the beginning each the year/month:
    if numel(dateRef) >= 3
        varargout{3}(ii) = 1;
        daysCurr = daysCurr - dateRef(3) + 1;
    end
    
    %Find output year:
    while daysCurr < daysIn(ii,:)
        if blCal == 0 && eomday(varargout{1}(ii),2) == 29
            daysAddYr = 366;
        elseif blCal == 1 || (blCal == 0 && eomday(varargout{1}(ii),2) == 28)
            daysAddYr = 365;
        elseif blCal == 2
            daysAddYr = 360;
        else
            error('days_2_date:unkownCal','Unkown calednar type input.');
        end
        
        %Add days to sum and one year to output if cumulative sum less than
        %or equal to daysIn
        if daysCurr + daysAddYr <=  daysIn(ii,:)
        	daysCurr = daysCurr + daysAddYr;
            varargout{1}(ii) = varargout{1}(ii) + 1;
        else
            break
        end
    end

    %Find month:
    while daysCurr < daysIn(ii,:)

        %May iterate past December, in which case must set to January 
        if varargout{2}(ii) == 13
           varargout{2}(ii) = 1;
           varargout{1}(ii) = varargout{1}(ii) + 1;
        end
        
        if varargout{2}(ii) == 2 && blCal == 1 %Select if February and 365 day calendar selected.
            daysAddMnth = 28;
        elseif blCal == 2
            daysAddMnth = 30;
        else
            daysAddMnth = eomday(varargout{1}(ii),varargout{2}(ii));
        end

        if daysCurr + daysAddMnth <=  daysIn(ii,:)
            daysCurr = daysCurr + daysAddMnth;
            varargout{2}(ii) = varargout{2}(ii) + 1;
        else
            break
        end
    end

    
    if varargout{2}(ii) > 12 %If month went over 12, transform to 1 
        varargout{2}(ii) = varargout{2}(ii) - 12;
        varargout{1}(ii) = varargout{1}(ii) + 1;
    elseif varargout{2}(ii) == 0 %If month is 0, transform  1 (probably cannot occur)
        varargout{1}(ii) = varargout{1}(ii) - 1;
        varargout{2}(ii) = 12;
    end


    %Find day of month:
    varargout{3}(ii) = daysIn(ii,:) - daysCurr + 1;
    
    %Record day fraction if nargout = 4
    if numel(varargout(:)) > 3
        for jj = 4 : numel(varargout(:))
            switch jj
                case 4
                    divisor = 24;
                case 5 || 6
                    divisor = 60;
            end
            int = floor(varargout{jj-1}(ii));
            varargout{jj}(ii) = (varargout{jj-1}(ii) - int)*divisor;
            varargout{jj-1}(ii) = int;
        end
    end
end


%[varargout{1},varargout{2},varargout{3},varargout{4}]

%Change 0 hr to 24, 0 min to 60, etc.
if numel(varargout(:)) >= 3
    for jj = numel(varargout(:)) : -1 : 3
        ind0 = find(varargout{jj} == 0);
%         switch jj
%             case 3
%                 
%             case 4
%                 varargout{jj}(ind0) = 24;
%             case 5 || 6
%                 varargout{jj}(ind0) = 60;
%         end
        
        for kk = 1 : numel(ind0)
            indX = zeros(jj,1);
            for ll = 1 : jj
                if varargout{ll}(ind0(kk)) == 1
                    indX(ll) = 1;
                end
            end
            %[varargout{1}(ind0(kk)), varargout{2}(ind0(kk)), varargout{3}(ind0(kk)), varargout{4}(ind0(kk))]
            if indX(jj-1) == 1
                if indX(jj-2) == 1
                    if jj - 3 > 0 && indX(jj-3) == 1
                        if jj - 4 > 0 && indX(jj-4) == 1
                            error('days_2_date:recursiveLimit',['The '...
                                'recursive limit has been reached for '...
                                'fixing dates that were inscribed as 0.']);
                        end
                    else
                        switch jj - 1
                            case 5
                                varargout{3}(ind0(kk)) = varargout{3}(ind0(kk)) - 1;
                                varargout{4}(ind0(kk)) = 24;
                                varargout{5}(ind0(kk)) = 60;
                                varargout{6}(ind0(kk)) = 60;
                            case 4
                                varargout{2}(ind0(kk)) = varargout{2}(ind0(kk)) - 1;
                                varargout{3}(ind0(kk)) = eomday(varargout{1}(ind0(kk)),varargout{2}(ind0(kk)));
                                varargout{4}(ind0(kk)) = 24;
                                varargout{5}(ind0(kk)) = 60;
                            case 3
                                varargout{1}(ind0(kk)) = varargout{1}(ind0(kk)) - 1;
                                varargout{2}(ind0(kk)) = 12;
                                varargout{3}(ind0(kk)) = eomday(varargout{1}(ind0(kk)),varargout{2}(ind0(kk)));
                                varargout{4}(ind0(kk)) = 24;
                            case 2
                                varargout{1}(ind0(kk)) = varargout{1}(ind0(kk)) - 1;
                                varargout{2}(ind0(kk)) = 12;
                                varargout{3}(ind0(kk)) = eomday(varargout{1}(ind0(kk)),varargout{2}(ind0(kk)));
                        end 
                    end
                else
                    switch jj - 1
                        case 5
                            varargout{4}(ind0(kk)) = varargout{4}(ind0(kk)) - 1;
                            varargout{5}(ind0(kk)) = 60;
                            varargout{6}(ind0(kk)) = 60;
                        case 4
                            varargout{3}(ind0(kk)) = varargout{3}(ind0(kk)) - 1;
                            varargout{4}(ind0(kk)) = 24;
                            varargout{5}(ind0(kk)) = 60;
                        case 3
                            varargout{2}(ind0(kk)) = varargout{2}(ind0(kk)) - 1;
                            varargout{3}(ind0(kk)) = eomday(varargout{1}(ind0(kk)),varargout{2}(ind0(kk)));
                            varargout{4}(ind0(kk)) = 24;
                        case 2
                            varargout{1}(ind0(kk)) = varargout{1}(ind0(kk)) - 1;
                            varargout{2}(ind0(kk)) = 12;
                            varargout{3}(ind0(kk)) = eomday(varargout{1}(ind0), varargout{2}(ind0));
                        case 1
                            varargout{1}(ind0(kk)) = varargout{1}(ind0(kk)) - 1;
                    end
                end
            else
                varargout{jj-1}(ind0(kk)) = varargout{jj-1}(ind0(kk)) - 1;
                switch jj
                    case 3
                        if varargout{2}(ind0) == 0 %This shouldn't be necessary here
                            varargout{2}(ind0) = 12;
                            varargout{1}(ind0) = varargout{1}(ind0) - 1;
                        end
                        varargout{3}(ind0) = eomday(varargout{1}(ind0), varargout{2}(ind0));
                    case 4
                        varargout{4}(ind0) = 24;
                    case 5 || 6
                        varargout{jj}(ind0) = 60;
                end
            end
            %find()
        end
%             indRep = 1; cntr = 1;
%             while indRep == 1
%                 if all(varargout{cntr+1:jj}(ind0(kk)) - 1 == 0)
%                     switch cntr
%                         case 1
%                             varargout{1}(ind0(kk)) = varargout{1}(ind0(kk)) + 1;
%                         case 2
%                             varargout{1}(ind0(kk)) = varargout{1}(ind0(kk)) + 1;
%                         case 6
%                             varargout{5}(ind0) = 60;
%                         case 5
%                             varargout{4}(ind0(kk)) = 24;
%                         case 4
%                             varargout{3}(ind0(kk)) = eomday(varargout{1}(ind0(kk)), varargout{2}(ind0(kk)));
%                         case 3
%                             varargout{2}(ind0(kk)) =
%                             
%                     end
%                 else
%                     varargout{jj-1}(ind0(kk)) = varargout{jj-1}(ind0(kk)) - 1;
%                     indRep = 0;
%                 end
%                 if cntr == jj
%                    indRep = 0; 
%                 else
%                     cntr = cntr + 1;
%                 end
%             end
%         end
    end
end
    
%Round (sometimes numerical errors crop up)
for ii = 1 : numel(varargout)
    varargout{ii} = round2(varargout{ii}, 3);
end

%Combine into single output:
if nargout <= 1
    for kk = 2 : numel(dateRef)
        varargout{1} = [varargout{1}, varargout{kk}];
    end
end


end