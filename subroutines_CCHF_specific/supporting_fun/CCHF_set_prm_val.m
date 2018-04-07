function cfAll = CCHF_set_prm_val(cfAll, cfKnown)



%Check to ensure no repeated parameter names:
for ii = 1 : numel(cfAll(:,1))
    indTest = setdiff((1:numel(cfAll(:,1))), ii);
    
    for jj = 1 : numel(indTest)
        if any(strcmpi(cfAll{ii,1}, cfAll{indTest(jj),1}))
            error('CCHFSetPrmVal:samePrmNm',['Two model parameters have the name ' cfAll{ii,1} '. One of them must be renamed in the appropriate process paramterization.'])
        end
    end
end

%Assign known/defined values to parameters
if numel(cfKnown(:)) > 0
    if numel(cfAll(1,:)) == 3 %This occurs for all types of runs except calibration
        for ii = 1 : numel(cfKnown(:,1))
            indCurr = find(strcmpi(cfAll(:,1), cfKnown{ii,2}) & strcmpi(cfAll(:,3), cfKnown{ii,1}) == 1);

            if ~isempty(indCurr)
                cfAll{indCurr, 2} = cfKnown{ii,3};

                disp(['The parameter ' cfKnown{ii,2} ' in function ' cfKnown{ii,1} ' has been set to ' num2str(cfAll{indCurr, 2}) '.']);
            else
               warning('CCHfsetPrm:missingParamter', ['The parameter ' cfKnown{ii,2} ...
                   ' used in ' cfKnown{ii,1} ' is being assigned a ' ...
                   'constant value but not not appear to be used in the prescribed module set.']); 
            end
        end
    elseif  numel(cfAll(1,:)) == 2 %This occurs for all types of runs except calibration
        for ii = 1 : numel(cfKnown(:,1))
            indCurr = find(strcmpi(cfAll(:,1), cfKnown{ii,2}));

            if ~isempty(indCurr)
                cfAll{indCurr, 2} = cfKnown{ii,3};

                disp(['The parameter ' cfKnown{ii,2} ' in function ' cfKnown{ii,1} ' has been set to ' num2str(cfAll{indCurr, 2}) '.']);
            else
               warning('CCHfsetPrm:missingParamter', ['The parameter ' cfKnown{ii,2} ...
                   ' used in ' cfKnown{ii,1} ' is being assigned a ' ...
                   'constant value but not not appear to be used in the prescribed module set.']); 
            end
        end
    else
       for ii = 1 : numel(cfKnown(:,1))
            indCurr = find(strcmpi(cfAll(:,1), cfKnown{ii,2}) & strcmpi(cfAll(:,5), cfKnown{ii,1}) == 1);

            if ~isempty(indCurr)
                cfAll{indCurr, 2} = cfKnown{ii,3};
                cfAll{indCurr, 3} = cfKnown{ii,3};
                cfAll{indCurr, 4} = cfKnown{ii,3};

                disp(['The parameter ' cfKnown{ii,2} ' in function ' cfKnown{ii,1} ' has been set to ' num2str(cfAll{indCurr, 2}) '.']);
            else
               warning('CCHfsetPrm:missingParamter', ['The parameter ' cfKnown{ii,2} ...
                   ' in function ' cfKnown{ii,1} ' appears to not be present in the prescribed module set.']); 
            end
        end
    end
end