function dt = time2sec(n,unit,varargin)
%Calculate number of seconds in given unit of time.

if regexpbl(unit,'month')
    if ~isempty(varargin(:))
        dt = 86400*eomday(varargin{1}(1),varargin{1}(2));
    else
        error('time2sec:curr_Date','The current year and month must be supplied since monthly time step used.');
    end
elseif regexpbl(unit,{'day','daily'})
    dt = 86400;
elseif regexpbl(unit,'hour')
    dt = 3600;
else
    error('time2sec:dt',[char(39) unit char(39) ' is an unknown time step.']);
end

dt = n*dt;