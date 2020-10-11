%
%noNameFunc   Post dist
%    
%

function [PM] = noNameFunc(PM, response)

    trial = length(PM.x);
    
    %Update PDF space: PDF = PostNTS(stimuli) 
    if response == 1
        PM.pdf = PM.PosteriorNextTrailSuccess(:,:,:,:,find(PM.stimRange == PM.xCurrent));
    else
        PM.pdf = PM.PosteriorNextTrialFailure(:,:,:,:,find(PM.stimRange == PM.xCurrent));
    end
    
    PM.pdf = PM.pdf./sum(sum(sum(sum(PM.pdf)))); % size[201 201 1 1]
    
    [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
    [minEntropy, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
    [trash, I] = PAL_findMax(PM.pdf);   
    PM.xCurrent = PM.stimRange(I);
    PM.x(trial+1) = PM.xCurrent;
    
    PM.threshold(trial) = sum(sum(sum(sum(PM.priorAlphas.*PM.pdf))));
    PM.slope(trial) = sum(sum(sum(sum(PM.priorBetas.*PM.pdf))));
    PM.guess(trial) = sum(sum(sum(sum(PM.priorGammas.*PM.pdf))));
    PM.lapse(trial) = sum(sum(sum(sum(PM.priorLambdas.*PM.pdf))));

  
    
   


