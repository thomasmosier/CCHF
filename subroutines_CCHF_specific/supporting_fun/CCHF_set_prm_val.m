function cfAll = CCHF_set_prm_val(cfAll, cfKnown)


if ~iscell(cfAll{1})
    blNotCell = 1;
    cfAll = {cfAll};
else
    blNotCell = 0;
end

%Check to ensure no repeated parameter names:
for mm = 1 : numel(cfAll)
    for ii = 1 : numel(cfAll{mm}(:,1))
        indTest = setdiff((1:numel(cfAll{mm}(:,1))), ii);

        for jj = 1 : numel(indTest)
            if any(strcmpi(cfAll{mm}{ii,1}, cfAll{mm}{indTest(jj),1}))
                error('CCHFSetPrmVal:samePrmNm',['Two model parameters have the name ' cfAll{mm}{ii,1} '. One of them must be renamed in the appropriate process paramterization.'])
            end
        end
    end

    %Assign known/defined values to parameters
    if numel(cfKnown(:)) > 0
        if numel(cfAll{mm}(1,:)) == 3 %This occurs for all types of runs except calibration
            for ii = 1 : numel(cfKnown(:,1))
                indCurr = find(strcmpi(cfAll{mm}(:,1), cfKnown{ii,2}) & strcmpi(cfAll{mm}(:,3), cfKnown{ii,1}) == 1);

                if ~isempty(indCurr)
                    cfAll{mm}{indCurr, 2} = cfKnown{ii,3};
                    
                    if mm == 1
                        disp(['The parameter ' cfKnown{ii,2} ' in subroutine ' cfKnown{ii,1} ' has been set to ' num2str(cfAll{mm}{indCurr, 2}) '.']);
                    end
%                 else
%                    warning('CCHfsetPrm:missingParamter', ['The parameter ' cfKnown{ii,2} ...
%                        ' used in ' cfKnown{ii,1} ' is being assigned a ' ...
%                        'constant value but not not appear to be used in the prescribed module set.']); 
                end
            end
        elseif  numel(cfAll{mm}(1,:)) == 2 %This occurs for all types of runs except calibration
            for ii = 1 : numel(cfKnown(:,1))
                indCurr = find(strcmpi(cfAll{mm}(:,1), cfKnown{ii,2}));

                if ~isempty(indCurr)
                    cfAll{mm}{indCurr, 2} = cfKnown{ii,3};

                    if mm == 1
                        disp(['The parameter ' cfKnown{ii,2} ' in subroutine ' cfKnown{ii,1} ' has been set to ' num2str(cfAll{mm}{indCurr, 2}) '.']);
                    end
%                 else
%                    warning('CCHfsetPrm:missingParamter', ['The parameter ' cfKnown{ii,2} ...
%                        ' used in ' cfKnown{ii,1} ' is being assigned a ' ...
%                        'constant value but not not appear to be used in the prescribed module set.']); 
                end
            end
        else
           for ii = 1 : numel(cfKnown(:,1))
                indCurr = find(strcmpi(cfAll{mm}(:,1), cfKnown{ii,2}) & strcmpi(cfAll{mm}(:,5), cfKnown{ii,1}) == 1);

                if ~isempty(indCurr)
                    cfAll{mm}{indCurr, 2} = cfKnown{ii,3};
                    cfAll{mm}{indCurr, 3} = cfKnown{ii,3};
                    cfAll{mm}{indCurr, 4} = cfKnown{ii,3};

                    if mm == 1
                        disp(['The parameter ' cfKnown{ii,2} ' in subroutine ' cfKnown{ii,1} ' has been set to ' num2str(cfAll{mm}{indCurr, 2}) '.']);
                    end
%                 else
%                    warning('CCHfsetPrm:missingParamter', ['The parameter ' cfKnown{ii,2} ...
%                        ' in function ' cfKnown{ii,1} ' appears to not be present in the prescribed module set.']); 
                end
            end
        end
    end
end

if blNotCell == 1
    if numel(cfAll(:)) == 1
        cfAll = cfAll{1};
    else
        error('CchfSetPrmVal:noCellNot1', ['The coef array has ' ...
            num2str(numel(cfAll(:))) ' parameter sets. But this ' ...
            'format option is only expected for arrays with one parameter set.']);
    end
end