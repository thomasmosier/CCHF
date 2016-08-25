function dt = time2hr(n,unit,varargin)
%Calculate number of seconds in given unit of time.

if regexpbl(unit,'month')
    if ~isempty(varargin(:))
        dt = 24*eomday(varargin{1}(1),varargin{1}(2));
    else
        error('time2hr:curr_Date','The current year and month must be supplied since monthly time step used.');
    end
elseif regexpbl(unit,{'day','daily'})
    dt = 24;
elseif regexpbl(unit,'hour')
    dt = 1;
else
    error('time2hr:dt',[char(39) unit char(39) ' is an unknown time step.']);
end

dt = n*dt;