%
%noNameFunc   Post dist
%    
%

function [PM] = UpdateFunc(PM, response)

    trial = length(PM.x);
    
% Step 7 
    % "Keep the posterior probability distribution from step 2 that 
    % corresponds to the completed trial" [Kontsevich]
    
    % pt+1(λ) = pt(λ|xt+1, rt+1)
    if response == 1
        PM.pdf = PM.PosteriorNextTrailSuccess(:,:,:,:,find(PM.stimRange == PM.xCurrent));
    else
        PM.pdf = PM.PosteriorNextTrialFailure(:,:,:,:,find(PM.stimRange == PM.xCurrent));
    end

% Step 8 
    % "Find a new estimate of the psychometric function based on the new 
    % posterior probability distribution pt+1(λ). The expected value of l provides the
    % answer:
    
    % λt+1 = sum of λ( λ*pt+1(λ) )
    PM.pdf = PM.pdf./sum(sum(sum(sum(PM.pdf)))); % size[201 201 1 1]

% Step 9
    [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
    [~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
    
    PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
    PM.x(trial+1) = PM.xCurrent;
    
    PM.threshold(trial) = sum(sum(sum(sum(PM.priorAlphas.*PM.pdf))));
    PM.slope(trial) = sum(sum(sum(sum(PM.priorBetas.*PM.pdf))));
    PM.guess(trial) = sum(sum(sum(sum(PM.priorGammas.*PM.pdf))));
    PM.lapse(trial) = sum(sum(sum(sum(PM.priorLambdas.*PM.pdf))));

  
    
   


