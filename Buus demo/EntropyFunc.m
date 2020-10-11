%
%EntropyFunc   Post dist
%    
%

function [minEntropy, newIntensityIndexPosition] = EntropyFunc(PosteriorNextTrailSuccess,PosteriorNextTrialFailure,pSuccessGivenx)

%TRIN 3: 
    %Ht(x,r)
    %success
    EntropyS = PosteriorNextTrailSuccess.*log(PosteriorNextTrailSuccess);
    EntropyS(isnan(EntropyS)) = 0;          %effectively defines 0.*log(0) to equal 0.
    for d = 1:4
        EntropyS = sum(EntropyS,d);    
    end
    EntropyS = -EntropyS;
    
    %Failure
    EntropyF = PosteriorNextTrialFailure.*log(PosteriorNextTrialFailure);
    EntropyF(isnan(EntropyF)) = 0;          %effectively defines 0.*log(0) to equal 0.
    for d = 1:4
        EntropyF = sum(EntropyF,d);    
    end
    EntropyF = -EntropyF;
    
%TRIN 4:  
    %E[Ht(x)]:      
    Entropy = EntropyS.*pSuccessGivenx + EntropyF.*(1-pSuccessGivenx); % [1,1,1,1,21]
    
%TRIN 5: 
    %xt+1 = arg min E[Ht(x)],
    [minEntropy, newIntensityIndexPosition] = min(squeeze(Entropy)); % minEntropy: the value, i; index position. 
    

