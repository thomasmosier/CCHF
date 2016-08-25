function nFit = fitness(obs, mod, strEval, varargin)


%Score for fitness metric evaluation that returns 0. 
nanPen = 10^3;

if numel(size(mod)) == 2
    %Ensure observed and modelled vectors have same orientation (probably
    %doesn't matter):
    if numel(obs(1,:)) == 1 && numel(obs(:,1)) ~= 1
       obs = obs'; 
    end
    if numel(mod(1,:)) == 1 && numel(mod(:,1)) ~= 1
       mod = mod'; 
    end
    
    %Remove values based on nan's in observed vector:
    mod(isnan(obs)) = [];
    obs(isnan(obs)) = [];

    %If optional input, remove values based on nan's in modeled vector:
    if ~isempty(varargin{:}) && regexpbl(varargin{1},'rmMod')
        obs(isnan(mod)) = [];
        mod(isnan(mod)) = [];
    end
end

if ischar(strEval)
   strEval = {strEval}; 
end

nFit = numel(strEval(:));
for mm = 1 : numel(strEval(:))
    
%     %Can have two evaluation methods, one = MODIS (only for 3D arrays) and
%     %second for all other cases:
%     if regexpbl(strEval{mm},'Parajka')
%         if numel(size(mod)) == 3 
%             strEval{mm} = 'Parajka';
%         elseif numel(size(mod)) == 2
%             indMOD = regexpi(strEval{mm},'Parajka');
%             strEval{mm}(indMOD:indMOD+4) = []; 
%         end
%     end

    if ~isequal(size(obs),size(mod))
       warning('fitness:unequalSize',['The model output and observations '...
           'to be compared have different sizes and are therefore not bieng compared.']); 
       nFit(mm) = nan;
       return
    end



    % %Make obserseved the same size as modelled so subtraction/add work.
    % if ~isequal(size(obs),size(mod))
    %     obs = obs(ones(size(mod,1),1),:);
    % end



    if regexpbl(strEval{mm},'bias') %BIAS:
        nFit(mm) = nanmean(obs(:) - mod(:));
    elseif regexpbl(strEval{mm},'RMSE') %ROOT MEAN SQUARE ERROR:
        nFit(mm) = sqrt(nanmean((mod(:) - obs(:)).^2));
    elseif regexpbl(strEval{mm},'abs') %MINIMIZE ABSOLUTE DIFFERENCE:
        nFit(mm) = sqrt(nansum((obs(:) - mod(:)).^2 , 2));
    elseif regexpbl(strEval{mm},{'MAE','mean absolute error'})
        nFit(mm) = nanmean( abs(mod(:) - obs(:)));
    elseif strcmpi(strEval{mm},'mape') %MINIMIZE MEAN ABSOLUTE PERCENT ERROR:
        indRem = find(obs == 0);
        indRem = union(indRem, find(isnan(obs)));
        indRem = union(indRem, find(isnan(mod)));
        obs(indRem) = [];
        mod(indRem) = [];

        nFit(mm) = nanmean( abs((mod(:) - obs(:))./obs(:)));
    elseif strcmpi(strEval{mm},'wmape') %MINIMIZE WEIGHTED MEAN ABSOLUTE PERCENT ERROR:
        indRem = find(obs == 0);
        indRem = union(indRem, find(isnan(obs)));
        indRem = union(indRem, find(isnan(mod)));
        obs(indRem) = [];
        mod(indRem) = [];
        nFit(mm) = nanmean( abs(obs(:)).*abs((mod(:) - obs(:))./obs(:)));
    elseif regexpbl(strEval{mm},'NSE') || regexpbl(strEval{mm},{'nash','sutcliffe'},'and')%MINIMIZE USING NASH-SUTCLIFFE EFFICIENCY:
        if numel(size(mod)) == 3
            %Calculate the KGE score over time at each grid location, then average
            %over domain
            nFitG = nan(size(squeeze(mod(1,:,:))));
            for jj = 1 : numel(mod(1,:,1))
                for ii = 1 : numel(mod(1,1,:))
                    indUse = intersect(find(~isnan(squeeze(obs(:,jj,ii)))), find(~isnan(squeeze(mod(:,jj,ii)))));

                    if ~isempty(indUse)
                        denomNse = nansum( (squeeze(obs(indUse,jj,ii)) - nanmean(squeeze(obs(indUse,jj,ii)),2).*ones(numel(squeeze(obs(indUse,jj,ii))),1)).^2);
                        nFitG(jj,ii) = nansum((squeeze(obs(indUse,jj,ii)) - squeeze(mod(indUse,jj,ii))).^2) ...
                            ./ denomNse;
                        if denomNse == 0
                            nFitG(jj,ii) = nan;
                        end 
                    else
                        nFitG(jj,ii) = nan;
                    end
                end
            end

            nFit(mm) = mean2d(nFitG);
        elseif numel(size(mod)) == 2
            denomNse = nansum( (obs - nanmean(obs,2)*ones(1, numel(obs(1,:)))).^2 , 2);
            nFit(mm) = nansum((obs - mod).^2, 2) ...
                ./ denomNse;
            if denomNse == 0
                nFit(mm) = nan;
            end 
        end
    elseif regexpbl(strEval{mm},'KGEr')
        if numel(size(mod)) == 3
            %Calculate the KGE score over time at each grid location, then average
            %over domain

            nFitG = nan(size(squeeze(mod(1,:,:))));
            for jj = 1 : numel(mod(1,:,1))
                for ii = 1 : numel(mod(1,1,:))
                    indUse = intersect(find(~isnan(squeeze(obs(:,jj,ii)))), find(~isnan(squeeze(mod(:,jj,ii)))));

                    if ~isempty(indUse)
                        covVar = cov(squeeze(obs(indUse,jj,ii)), squeeze(mod(indUse,jj,ii)));

                        nFitG(jj,ii) = abs(covVar(2)/(std(squeeze(obs(indUse,jj,ii)))*std(squeeze(mod(indUse,jj,ii))))-1);
                    else
                        nFitG(jj,ii) = nan;
                    end
                end
            end

            nFit(mm) = mean2d(nFitG);
        elseif numel(size(mod)) == 2
            covVar = cov(obs,mod);
            nFit(mm) = abs(covVar(2)/(std(obs)*std(mod))-1);
        end
    elseif regexpbl(strEval{mm},'KGEs')
        if numel(size(mod)) == 3
            %Calculate the KGE score over time at each grid location, then average
            %over domain

            nFitG = nan(size(squeeze(mod(1,:,:))));
            for jj = 1 : numel(mod(1,:,1))
                for ii = 1 : numel(mod(1,1,:))
                    indUse = intersect(find(~isnan(squeeze(obs(:,jj,ii)))), find(~isnan(squeeze(mod(:,jj,ii)))));

                    if ~isempty(indUse)
                        nFitG(jj,ii) = abs(std(mod(indUse,jj,ii))/std(obs(indUse,jj,ii))-1);
                    else
                        nFitG(jj,ii) = nan;
                    end
                end
            end

            nFit(mm) = mean2d(nFitG);
        elseif numel(size(mod)) == 2
            nFit(mm) = abs(std(mod)/std(obs)-1);
        end
    elseif regexpbl(strEval{mm},'KGEb')
        if numel(size(mod)) == 3
            %Calculate the KGE score over time at each grid location, then average
            %over domain

            nFitG = nan(size(squeeze(mod(1,:,:))));
            for jj = 1 : numel(mod(1,:,1))
                for ii = 1 : numel(mod(1,1,:))
                    indUse = intersect(find(~isnan(squeeze(obs(:,jj,ii)))), find(~isnan(squeeze(mod(:,jj,ii)))));

                    if ~isempty(indUse)
                        nFitG(jj,ii) = abs(nanmean(mod(indUse,jj,ii))/nanmean(obs(indUse,jj,ii))-1);
                    else
                        nFitG(jj,ii) = nan;
                    end
                end
            end

            nFit(mm) = mean2d(nFitG);
        elseif numel(size(mod)) == 2
            nFit(mm) = abs(nanmean(mod)/nanmean(obs)-1);
        end
    elseif (regexpbl(strEval{mm},'KGE') || regexpbl(strEval{mm},{'kling','gupta'},'and')) %&& ~regexpbl(strEval{mm}, '3D')
        if numel(size(mod)) == 3
            %Calculate the KGE score over time at each grid location, then average
            %over domain

            nFitG = nan(size(squeeze(mod(1,:,:))));
            for jj = 1 : numel(mod(1,:,1))
                for ii = 1 : numel(mod(1,1,:))
                    indUse = intersect(find(~isnan(squeeze(obs(:,jj,ii)))), find(~isnan(squeeze(mod(:,jj,ii)))));

                    if ~isempty(indUse)
                        %disp([num2str(jj) ', ' num2str(ii) ': ' num2str(numel(indUse))])

                        %[r,p] = corrcoef(squeeze(obs(indUse,jj,ii)), squeeze(mod(indUse,jj,ii)));

                        covVar = cov(squeeze(obs(indUse,jj,ii)), squeeze(mod(indUse,jj,ii)));
                        %See definition of the KGE and linear correlation coefficient (can be decomposed into covariance and SD)
                    %     nFit(mm) = sqrt((corrcoef(2)-1)^2 + (std(mod)/std(obs)-1)^2 + (nanmean(mod)/nanmean(obs)-1)^2);

                        nFitG(jj,ii) = sqrt((covVar(2)/(std(squeeze(obs(indUse,jj,ii)))*std(squeeze(mod(indUse,jj,ii))))-1)^2 ...
                            + (std(squeeze(mod(indUse,jj,ii)))/std(squeeze(obs(indUse,jj,ii)))-1)^2 ...
                            + (nanmean(squeeze(mod(indUse,jj,ii)))/nanmean(squeeze(obs(indUse,jj,ii)))-1)^2);
                    else
                        nFitG(jj,ii) = nan;
                    end
                end
            end

            nFit(mm) = mean2d(nFitG);
        elseif numel(size(mod)) == 2
            covVar = cov(obs,mod);
            %See definition of the KGE and linear correlation coefficient (can be decomposed into covariance and SD)
        %     nFit(mm) = sqrt((corrcoef(2)-1)^2 + (std(mod)/std(obs)-1)^2 + (nanmean(mod)/nanmean(obs)-1)^2);
            nFit(mm) = sqrt((covVar(2)/(std(obs)*std(mod))-1)^2 + (std(mod)/std(obs)-1)^2 + (nanmean(mod)/nanmean(obs)-1)^2);
        end
    elseif regexpbl(strEval{mm},'Pearson') || regexpbl(strEval{mm},{'Pearson','correlation','coefficient'},'and')
        if numel(size(mod)) == 3

            nFitG = nan(size(squeeze(mod(1,:,:))));
            for jj = 1 : numel(mod(1,:,1))
                for ii = 1 : numel(mod(1,1,:))
                    indUse = intersect(find(~isnan(squeeze(obs(:,jj,ii)))), find(~isnan(squeeze(mod(:,jj,ii)))));

                    if ~isempty(indUse)
                        [r,p] = corrcoef(obs(indUse,jj,ii), mod(indUse,jj,ii));
                        nFitG(jj,ii) = [r(2), p(2)];
                    else
                        nFitG(jj,ii) = nan;
                    end
                end
            end

            nFit(mm) = mean2d(nFitG);
        elseif numel(size(mod)) == 2
            [r,p] = corrcoef(obs,mod);
            nFit(mm) = [r(2), p(2)];
        end
    elseif regexpbl(strEval{mm},{'Parajka','MODIS'})
        %Find number of pixels that are glaciated or have permanent snow cover
        %(i.e. always 100 or nan in MODIS observations)
        subGlac = zeros(0,2);
        nGlacThresh = 0.01*numel(obs(:,1,1));
        for jj = 1 : numel(obs(1,:,1))
            for ii = 1 : numel(obs(1,1,:))
                if sum2d(~isnan(obs(:,jj,ii)) & obs(:,jj,ii) ~= 100) < nGlacThresh
                    subGlac(end+1,:) = [jj,ii];

                    mod(:,jj,ii) = nan;
                    obs(:,jj,ii) = nan;
                end
            end
        end

        %Parajka, J., & Blöschl, G. (2008). The value of MODIS snow cover data 
        %in validating and calibrating conceptual hydrologic models. Journal of 
        %Hydrology, 358(3?4), 240-258. http://doi.org/10.1016/j.jhydrol.2008.06.006
        
        %Find number of images that are not all nan:
        mCntr = 0;
        for kk = 1 : numel(obs(:,1,1))
            if ~all2d(isnan(squeeze(obs(kk,:,:))))
               mCntr = mCntr + 1; 
            end
        end
        
        lPix = numel(obs(1,:,:)) - numel(subGlac(:,1));

        thrshSWE = 10/100; %units = m
        thrshSCA = 10; %units = percent
        wghtOvr = 5; %Weighting factor for Over
        wghtUnd = 5; %Weighting factor for under
        
        EOvr = nan(size(squeeze(obs(1,:,:))));
        EUnd = nan(size(squeeze(obs(1,:,:))));
        for jj = 1 : numel(obs(1,:,1))
            for ii = 1 : numel(obs(1,1,:))
                if ~ismember([jj,ii],subGlac,'rows')
                    indNNan = find(~isnan(mod(:,jj,ii)) & ~isnan(obs(:,jj,ii)));

                    indYObs = find(obs(indNNan,jj,ii) > thrshSCA);
                    indNObs = find(obs(indNNan,jj,ii) == 0);
                    indYMod = find(mod(indNNan,jj,ii) > thrshSWE);
                    indNMod = find(mod(indNNan,jj,ii) == 0);

                    EOvr(jj,ii) = numel(intersect(indYMod, indNObs));
                    EUnd(jj,ii) = numel(intersect(indNMod, indYObs));
                end
            end
        end
        EOvrAvg = sum2d(EOvr)/(mCntr*lPix);
        EUndAvg = sum2d(EUnd)/(mCntr*lPix);

        nFit(mm) = wghtOvr*EOvrAvg + wghtUnd*EUndAvg;

    else
        error('Unkown fitness ranking method specified.');
    end

    if isempty(nFit(mm) )
       nFit(mm) = nan; 
    end

    if isnan(nFit(mm) )
        nFit(mm) = nanPen;
    end
end