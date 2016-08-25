function [transectData, hdr] = climate_crossection(climateStr, transectDd)
crop = 'out';   %Out is preferable for longitude, but maybe not for latitude

if regexpi(climateStr,'.asc')
    [transectData, hdr] = read_ESRI(climateStr);
elseif regexpi(climateStr,'.bil')
    %Read and load header of WorldClim binary file:
    [pathClimate, nameClimate, extClimate] = fileparts(climateStr); 
    hdr = read_wc_bin_hdr(fullfile(pathClimate, [nameClimate '.hdr']));    
end

transectDd = [transectDd(1), transectDd(2), (transectDd(3)-hdr(5)),(transectDd(3)+hdr(5))];

if regexpi(climateStr,'.bil')
    %Load WC binary file:
    [transectData, hdr, ~, ~] = read_wc_bin_data(climateStr, hdr, transectDd, crop);  
elseif regexpi(climateStr,'.asc')
    [bndsDd, bndsInd] = adjust_bounds(hdr, transectDd, crop);
    
    if bndsInd(1) > hdr(1) || bndsInd(2) > hdr(1)
        error('The transect longitudinal bounds extend beyond the geographic scope of the file chosen.');
    elseif bndsInd(3) > hdr(2) || bndsInd(4) > hdr(2)
        error(['The transect' char(39) 's latitude is not within the bounds of one of the files chosen.']);
    end

    transectData = transectData(bndsInd(4) : bndsInd(3) , bndsInd(1) : bndsInd(2));
    transectData(transectData == hdr(6)) = NaN;
    
    hdr(1) = length(transectData(1,:));
    hdr(2) = length(transectData(:,1));
    hdr(3) = bndsDd(1) - hdr(5)/2;
    hdr(4) = bndsDd(3) - hdr(5)/2;
end


end
