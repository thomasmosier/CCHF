function coefOut = CCHF_ld_prms(pathPrm, sMeta)

if regexpbl(pathPrm,'.txt')
    [cfTemp, ~] = read_CCHF_coef(pathPrm);

    %Check that loaded coefficients are same as those needed for current
    %model:
    if ~regexpbl(cfTemp(:,1), sMeta.coefAtt(:,1),'and')
        %Find missing parameter:
        cfBl = zeros(numel(sMeta.coefAtt(:,1)),1);
        for ii = 1 : numel(sMeta.coefAtt(:,1))
            if regexpbl(cfTemp(:,1),sMeta.coefAtt{ii,1})
                cfBl(ii) = 1;
            end
        end
        indCf = find(cfBl == 0);
        if ~isempty(indCf)
            nmQuit = -9999;
            uiwait(msgbox(sprintf(['The loaded parameter set does not '...
                'include ' num2str(numel(indCf)) ' parameters needed '...
                'for the ' char(39) sMeta.strModule char(39) ' model '...
                'version.' char(10) 'You will be prompted to enter '...
                'values for these parameters. If you wish to cancel '...
                'this model run, enter ' num2str(nmQuit) ' for one of the values.']),...
                '(Click OK to Proceed)','modal'));

            for ii = 1 : numel(indCf)
               cfCurr = input(['Enter a value for ' char(39) ...
                   sMeta.coefAtt{indCf(ii),1} char(39) '. The bounds are ' ...
                   num2str(sMeta.coefAtt{indCf(ii),2}) ' thru ' ...
                   num2str(sMeta.coefAtt{indCf(ii),3}) ' and the parameter is '...
                   'used in ' char(39) sMeta.coefAtt{indCf(ii),5} char(39) ':' char(10)]); 
               cfTemp(end+1,:) = {sMeta.coefAtt{indCf(ii),1}, cfCurr};
               if cfCurr == nmQuit
                   error('CCHF_main:missingParam',['The exit command '...
                       'was entered during the process of entinering a'...
                       ' missing parameter value.']);
               end
            end
            clear ii
        end
    end

    coefOut = cfTemp;
elseif regexpbl(pathPrm,'.csv')
    prmArray = csvread(pathPrm,1,0);
    fidPrm = fopen(pathPrm,'r');
    varPrm = textscan(fidPrm,'%s',numel(prmArray(1,:)),'Delimiter',',');  
        varPrm = varPrm{1};
    fclose(fidPrm);

    if numel(varPrm(:)) - 1 == numel(sMeta.coefAtt(:,1))
        varPrm = varPrm(1:end-1)'; %Last column is evaluation metric
        prmArray = prmArray(:,1:end-1);
    elseif numel(varPrm(:)) ~= numel(sMeta.coefAtt(:,1))
        error('CCHF_main:nParam',['The number of parameters loaded '...
            'is not equal to the number required by this model formulation.']);
    end

    %Check that all parameters present:
    if ~regexpbl(varPrm, sMeta.coefAtt(:,1))
        error('CCHF_main:missingParam',['There is a missing '...
            'parameter in the current list of parameter sets.']);
    end

    %Check that order of parameters is correct:
    for ii = 1 : numel(varPrm)
        if ~regexpbl(varPrm{ii}, sMeta.coefAtt{ii,1})
            error('CCHF_main:outoforderParam',[varPrm{ii} ' in the '...
                'loaded parameter set is out of order relative to '...
                'the parameter set identified for the current model '...
                'formulation.']);
        end
    end
    clear ii
    
    %Populate output cell array, where each entry is a cell array of
    %parameter names and values
    szPrmIn = size(prmArray);
    coefOut = cell([szPrmIn(1), 1]);
    [coefOut{:}] = deal(cell([szPrmIn(2), 2]));
    
    for mm = 1 : numel(coefOut(:))
        coefOut{mm}(:,1) = varPrm(:);
        coefOut{mm}(:,2) = num2cell(prmArray(mm,:))';
    end
else
    [~, ~, extCf] = fileparts(pathPrm);
    error('CCHF_main:unknownCoefFormat',['The coeffienct file is in a ' ...
        char(39) extCf char(39) ' format, which has not been programmed for.']);
end