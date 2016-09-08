function fileDiary = downscale_diary(dirOutput,fileAppend)
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


%Open diary 
%If 'fileAppend' does not exist or is empty, use filename based upon 
%current time.
%If 'fileAppend' is not empty, use that file for creating the diary.

if nargin == 1 || isempty(fileAppend)
    tDiary = clock;
    fileDiary = fullfile(dirOutput,['processing_log_' ...
        num2str(tDiary(5)) '-' num2str(tDiary(4)) '-' ...
        num2str(tDiary(3)) '-' num2str(tDiary(2)) '-' ...
        num2str(tDiary(1)) '.txt']);
else
    if ~isempty(regexpi(dirOutput,fileAppend))
        fileDiary = fileAppend;
    else
        [~,file,ext] = fileparts(fileAppend);
        fileDiary = fullfile(dirOutput,[file,ext]);
    end

end

diary(fileDiary);