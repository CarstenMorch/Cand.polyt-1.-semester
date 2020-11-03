%
%BruteForceSearchFunc   Find maximum of of 4D parameters grid.
%
%   syntax: [parameters ] = ...
%       PAL_PFML_BruteForceFit(StimLevels, NumPos, OutOfNum, ...
%       searchGrid, PF)
%
% Gives sames results as PAL_PFML_BruteForceFit
%
%Latest edited: 
%   06-10-2020


function [ parameters ] = BruteForceSearchFunc(stimLevels,NumPos, OutOfNum, searchGrid, PF)
    % Create a likelihood grid 
    likelihoodGrid = zeros(length(searchGrid.alpha),length(searchGrid.beta),length(searchGrid.gamma),length(searchGrid.lambda));
    
    % Find likelihood grid 
    for a = 1:length(searchGrid.alpha)
        for b = 1:length(searchGrid.beta)
            for g = 1:length(searchGrid.gamma)
                for L = 1:length(searchGrid.lambda) 
                    for sLevel = 1:length(stimLevels)
                        temp(sLevel) = (PF([searchGrid.alpha(a), searchGrid.beta(b), searchGrid.gamma(g), searchGrid.lambda(L)], stimLevels(sLevel))^NumPos(sLevel))*((1-PF([searchGrid.alpha(a), searchGrid.beta(b), searchGrid.gamma(g), searchGrid.lambda(L)], stimLevels(sLevel)))^(OutOfNum(sLevel)-NumPos(sLevel)));
                    end
                    likelihoodGrid(a,b,g,L)=prod(temp); 
                    clear temp; 
                end
            end
        end 
    end 
    
    % Find maximum value in likelihood grid 
    Amax = max(likelihoodGrid(:));                                           % Maximum Value
    Idx = find(likelihoodGrid(:) == Amax);                                   % Returns Linear Indices To All ‘Amax’ Values
    [r,c,dim3,dim4] = ind2sub(size(likelihoodGrid), Idx); 
    
    % Return found values 
    parameters = [searchGrid.alpha(r) searchGrid.beta(c) searchGrid.gamma(dim3) searchGrid.lambda(dim4)];
end 
