function coefAll = CCHF_set_prm_val(coefAll, knownCoef)

if numel(knownCoef(:)) > 0
    if numel(coefAll(1,:)) == 3 %This occurs for all types of runs except calibration
        for ii = 1 : numel(knownCoef(:,1))
            indCurr = find(strcmpi(coefAll(:,1), knownCoef{ii,2}) & strcmpi(coefAll(:,3), knownCoef{ii,1}) == 1);

            if ~isempty(indCurr)
                coefAll{indCurr, 2} = knownCoef{ii,3};

                disp(['The parameter ' knownCoef{ii,2} ' in function ' knownCoef{ii,1} ' has been set to ' num2str(coefAll{indCurr, 2}) '.']);
            else
               warning('CCHfsetPrm:missingParamter', ['The parameter ' knownCoef{ii,2} ...
                   ' in function ' knownCoef{ii,1} ' appears to not be present in the prescribed module set.']); 
            end
        end
    else
       for ii = 1 : numel(knownCoef(:,1))
            indCurr = find(strcmpi(coefAll(:,1), knownCoef{ii,2}) & strcmpi(coefAll(:,5), knownCoef{ii,1}) == 1);

            if ~isempty(indCurr)
                coefAll{indCurr, 2} = knownCoef{ii,3};
                coefAll{indCurr, 3} = knownCoef{ii,3};
                coefAll{indCurr, 4} = knownCoef{ii,3};

                disp(['The parameter ' knownCoef{ii,2} ' in function ' knownCoef{ii,1} ' has been set to ' num2str(coefAll{indCurr, 2}) '.']);
            else
               warning('CCHfsetPrm:missingParamter', ['The parameter ' knownCoef{ii,2} ...
                   ' in function ' knownCoef{ii,1} ' appears to not be present in the prescribed module set.']); 
            end
        end
    end
end