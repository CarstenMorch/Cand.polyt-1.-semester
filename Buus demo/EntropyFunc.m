%
%EntropyFunc   Post dist
%    
%

function [minEntropy, newIntensityIndexPosition] = EntropyFunc(PosteriorNextTrailSuccess,PosteriorNextTrialFailure,pSuccessGivenx)

% Step 3: 
    % "Estimate the entropy of the probability density function over the 
    % space of psychometric functions, given that at the next trial a test 
    % of intensity x will produce the response r" [Kontsevich]
    
    % Ht(x,r)
    % success
    EntropyS = PosteriorNextTrailSuccess.*log(PosteriorNextTrailSuccess);
    EntropyS(isnan(EntropyS)) = 0;          %effectively defines 0.*log(0) to equal 0.
    for d = 1:4
        EntropyS = sum(EntropyS,d);    
    end
    EntropyS = -EntropyS;
    
    % failure
    EntropyF = PosteriorNextTrialFailure.*log(PosteriorNextTrialFailure);
    EntropyF(isnan(EntropyF)) = 0;          %effectively defines 0.*log(0) to equal 0.
    for d = 1:4
        EntropyF = sum(EntropyF,d);    
    end
    EntropyF = -EntropyF;
    
% Step 4:  
    % "Estimate the expected entropy for each test intensity x" 
    % [Kontsevich]
    
    % E[Ht(x)]:      
    Entropy = EntropyS.*pSuccessGivenx + EntropyF.*(1-pSuccessGivenx); % [1,1,1,1,21]
    
% Step 5:
    % "Find the test intensity that has the minimum expected entropy"
    % [Kontsevich]
    
    % xt+1 = arg min E[Ht(x)]:
    [minEntropy, newIntensityIndexPosition] = min(squeeze(Entropy)); % minEntropy: the value, i; index position. 
    

