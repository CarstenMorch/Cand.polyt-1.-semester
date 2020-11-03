%
%SingleValueLikelihood   Find maximum of of 4D parameters grid.
%
%   syntax: [parameters ] = ...
%       PAL_PFML_BruteForceFit(StimLevels, NumPos, OutOfNum, ...
%       searchGrid, PF)
%

function [ funcValue ] = SingleValueLikelihoodFunc(stimLevels,NumPos, OutOfNum, Parameter, PF)

    for sLevel = 1:length(stimLevels)
        temp(sLevel) = (PF([Parameter(1),Parameter(2),Parameter(3),Parameter(4)], stimLevels(sLevel))^NumPos(sLevel))*((1-PF([Parameter(1),Parameter(2),Parameter(3),Parameter(4)], stimLevels(sLevel)))^(OutOfNum(sLevel)-NumPos(sLevel)));
    end
    funcValue  = prod(temp); 

end 
    
   