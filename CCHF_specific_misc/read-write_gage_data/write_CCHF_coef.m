function write_CCHF_coef(path,hdr,coef,varargin)

% if isnumeric(coef)
%     
% end
permission = 'wt';
% strHdr = blanks(0);
if ~isempty(varargin)
    for ii = 1 : numel(varargin)
        switch varargin{ii}
            case 'permissions'
                permission = varargin{ii+1};
%             case 'hdr'
%                 strHdr = varargin{ii+1};
        end
    end
end
        

    
sNum = '%.6g';

fid = fopen(path, permission);


%Write Standard Header:
if iscell(hdr)
    for ii = 1 : numel(hdr(:,1))
        if numel(hdr(1,:)) == 1
            fprintf(fid,'%s\n',hdr{ii,:});
        elseif numel(hdr(1,:)) == 2
            if isnumeric(hdr(1,2))
                fprintf(fid, ['%s\t' sNum '\n'],hdr{ii,:});
            else
                fprintf(fid,'%s\t%s\n',hdr{ii,:});
            end
        else
            error('write_CCHF_coef:hdrNumel',['The input ' char(39) ...
                'hdr' char(39) ' has ' num2str(numel(hdr(1,:))) ...
                ' columns.  This has not been programmed for.']);
        end
    end
elseif ischar(hdr)
    fprintf(fid,'%s\n',hdr);
else
    error('write_CCHF_coef:hdrType',['The input ' char(39) 'hdr' ...
        char(39) ' is in a format that has not been programmed for.']);
end


%Write parameters:
for jj = 1 : numel(coef(:,1))
    if numel(coef(1,:)) == 1
        if ischar(coef{1,1})
        	fprintf(fid, '%s\n',coef{jj,:});
        elseif isnumeric(coef{1,1})
            fprintf(fid, [Num '\n'],coef{jj,:});
        else
            error('write_CCHF_coef:coefType',['The input ' char(39) 'coef' ...
                char(39) ' is in a format that has not been programmed for.']);
        end
    elseif numel(coef(1,:)) == 2
        if ischar(coef{1,1}) && isnumeric(coef{1,2})
            fprintf(fid, ['%s\t' sNum '\n'],coef{jj,:});
        elseif isnumeric(coef{1,1}) && ischar(coef{1,2})
            fprintf(fid, [sNum '\t%s\n'],coef{jj,:});
        else
            error('write_CCHF_coef:coefType',['The input ' char(39) 'coef' ...
                char(39) ' is in a format that has not been programmed for.']);
        end
    elseif numel(coef(1,:)) == 5
        if ischar(coef{1,1}) && isnumeric(coef{1,2}) && isnumeric(coef{1,3}) && isnumeric(coef{1,4}) && ischar(coef{1,5})
            fprintf(fid, ['%s\t' sNum '\t' sNum '\t' sNum '\t%s\n'],coef{jj,:});
        else
            error('write_CCHF_coef:coefType',['The input ' char(39) 'coef' ...
                char(39) ' is in a format that has not been programmed for.']);
        end
    end
end


fclose(fid);