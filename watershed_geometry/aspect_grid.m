function aspect = aspect_grid(fdr)
%'aspect' is the slope's azimuth (degrees clockwise from North)



aspect = nan(size(fdr));

for jj = 2 : numel(fdr(:,1)) - 1
	for ii = 2 : numel(fdr(1,:)) - 1     
        ind = fdr_inv_mult(fdr(jj,ii));
        
        out = nan(numel(ind),1);

        for kk = 1 : numel(ind)
            switch ind(kk)
            %E  =   1;      8
            %SE =   2;      9
            %S  =   4;      6
            %SW =   8;      3
            %W  =  16;      2
            %NW =  32;      1
            %N  =  64;      4
            %NE = 128;      7
                case 8
                    out = 90;
                case 9
                    out = 135;
                case 6
                    out = 180;
                case 3
                    out = 225;
                case 2
                    out = 270;
                case 1
                    out = 315;
                case 4
                    out = 0;
                case 7
                    out = 45;
                case 0
                    out = 5;
            end
        end

        aspect(jj,ii) = nanmean(out);
	end
end