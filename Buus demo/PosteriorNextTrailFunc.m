%
%PosteriorNextTrailFunc   Post dist
%   copy of PAL_AMPM_PosteriorTplus1 
%

function [PosteriorNextTrailSuccess,PosteriorNextTrialFailure,pSuccessGivenx] = PosteriorNextTrailFunc(pdf,LUT)
    % obs: LUT = p(r|lambda,x)
   
% Step 1: 
    % "Calculate probability of getting response r after presenting test x 
    % at the next trial" [Kontsevich]
    
    % pt(λ)
    pdf5D = repmat(pdf, [1 1 1 1 size(LUT,5)]); % 201 201 --> 201 201 1 1 21 
    
    % pt(λ)p(r|λ,x):
    PosteriorNextTrailSuccess = pdf5D.*LUT; % size[201 201 1 1 21] 
    PosteriorNextTrialFailure = pdf5D-PosteriorNextTrailSuccess;
    
    % p(succes|x) = sum of λ ( p(succes|λ,x)*pt(λ) ): 
    DenominatorS = squeeze(sum(sum(sum(sum(PosteriorNextTrailSuccess,1),2),3),4)); %the probability of stimulus in relation to all parameters, return the 1,1,1,1,21 --> 21,1
    DenominatorS = repmat(DenominatorS, [1 size(pdf5D,1) size(pdf5D,2) size(pdf5D,3) size(pdf5D,4)]); %[21,1] --> [21, 201, 201]    
    DenominatorS = permute(DenominatorS, [2 3 4 5 1]); %[21, 201, 201] --> [201 201 1 1 21]
      
    DenominatorF = squeeze(sum(sum(sum(sum(PosteriorNextTrialFailure,1),2),3),4));
    DenominatorF = repmat(DenominatorF, [1 size(pdf5D,1) size(pdf5D,2) size(pdf5D,3) size(pdf5D,4)]);
    DenominatorF = permute(DenominatorF, [2 3 4 5 1]);
    
% Step 2: 
    % "Estimate by Bayes’ rule the posterior probability of each 
    % psychometric function given that the next trial will produce the 
    % response r to the test of the intensity x." [Kontsevich]
    
    % In other words: The probability of alpha and beta at intensity x 
    % divided by the probability of the summed probability of all alpha and
    % beta values at intensity x.
    
    % p_t (λ|x,r): 
    PosteriorNextTrailSuccess = PosteriorNextTrailSuccess./DenominatorS; 
    PosteriorNextTrialFailure = PosteriorNextTrialFailure./DenominatorF;
    
    % Needed in step 4
    % pt(success|x):
    pSuccessGivenx = sum(sum(sum(sum(pdf5D.*LUT,1),2),3),4); 


