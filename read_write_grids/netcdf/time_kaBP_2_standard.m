function gcmDateVec = time_kaBP_2_standard(gcmTime,gcmRef,cal)
%Function only necessary when NetCDF units are 'ka BP' (thousands of years
%before present).  Assumes monthly data and that calendar is 'no leap'

if regexpbl(cal,{'no','leap'},'and')
    gcmDateVec = nan(numel(gcmTime),3);
    gcmDateVec(:,1) = gcmRef(1) + ceil(10^3*gcmTime); 
    gcmDateVec(:,2) = rem((1:numel(gcmTime)),12);
    gcmDateVec(gcmDateVec(:,2) == 0, 2) = 12;
    gcmDateVec(:,3) = 15;
else
   error('time_kaBP_2_standard:calendar',['No case has been written for calendar of type ' cal]);
end


end