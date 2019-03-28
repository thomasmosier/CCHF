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

function t = rmfield_x(s,field)
%RMFIELD Remove fields from a structure array.
%   S = RMFIELD(S,'field') removes the specified field from the
%   m x n structure array S. The size of input S is preserved.
%
%   S = RMFIELD(S,FIELDS) removes more than one field at a time
%   when FIELDS is a character array or cell array of strings.  The
%   changed structure is returned. The size of input S is preserved.
%
%   See also SETFIELD, GETFIELD, ISFIELD, FIELDNAMES.

%   Copyright 1984-2008 The MathWorks, Inc.

%--------------------------------------------------------------------------------------------
% handle input arguments
if ~isa(s,'struct') 
    error(message('MATLAB:rmfield:Arg1NotStructArray')); 
end
if ~ischar(field) && ~iscellstr(field)
   error(message('MATLAB:rmfield:FieldnamesNotStrings'));
elseif ischar(field)
   field = cellstr(field); 
end

field = unique(field); %In case repeat entries

% get fieldnames of struct
f = fieldnames(s);

% Determine which fieldnames to delete.
idxremove = [];
for i=1:length(field)
  j = find(strcmp(field{i},f) == true);
  if isempty(j)
    if length(field{i}) > namelengthmax
      error(message('MATLAB:rmfield:FieldnameTooLong', field{ i }));
%     else
%       error(message('MATLAB:rmfield:InvalidFieldname', field{ i }));
    end
  end
  idxremove = [idxremove;j];
end

% set indices of fields to keep
idxkeep = 1:length(f);
idxkeep(idxremove) = [];

% remove the specified fieldnames from the list of fieldnames.
f(idxremove,:) = [];

% convert struct to cell array
c = struct2cell(s);

% find size of cell array
sizeofarray = size(c);
newsizeofarray = sizeofarray;

% adjust size for fields to be removed
newsizeofarray(1) = sizeofarray(1) - length(idxremove);

% rebuild struct
t = cell2struct(reshape(c(idxkeep,:),newsizeofarray),f);

%--------------------------------------------------------------------------------------------
