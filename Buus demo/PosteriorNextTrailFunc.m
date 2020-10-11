%
%PosteriorNextTrailFunc   Post dist
%   copy of PAL_AMPM_PosteriorTplus1 
%

function [PosteriorNextTrailSuccess,PosteriorNextTrialFailure,pSuccessGivenx] = PosteriorNextTrailFunc(pdf,LUT)
    % osb: LUT = p(r|lambda,x)
   
% TRIN 1: 
    pdf5D = repmat(pdf, [1 1 1 1 size(LUT,5)]); % 201 201 --> 201 201 1 1 21 
    
    % pt(lambda)p(r|lambda,x):
    PosteriorNextTrailSuccess = pdf5D.*LUT; % size[201 201 1 1 21] 
    PosteriorNextTrialFailure = pdf5D-PosteriorNextTrailSuccess;
    
    % p(succes|x) = sum( p(succes|lambda,x) * pt(lambda), Lambda: 
    DenominatorS = squeeze(sum(sum(sum(sum(PosteriorNextTrailSuccess,1),2),3),4)); %sandsynligheden for x ift all parameter, return the 1,1,1,1,21 --> 21,1
    DenominatorS = repmat(DenominatorS, [1 size(pdf5D,1) size(pdf5D,2) size(pdf5D,3) size(pdf5D,4)]); %[21,1] --> [21, 201, 201]    
    DenominatorS = permute(DenominatorS, [2 3 4 5 1]); %[21, 201, 201] --> [201 201 1 1 21]
    
    DenominatorF = squeeze(sum(sum(sum(sum(PosteriorNextTrialFailure,1),2),3),4));
    DenominatorF = repmat(DenominatorF, [1 size(pdf5D,1) size(pdf5D,2) size(pdf5D,3) size(pdf5D,4)]);
    DenominatorF = permute(DenominatorF, [2 3 4 5 1]);
    
% TRIN 2: 
    % Sandsynligheden for alpha og beta ved intensitet x divideret med 
    % sandsynligheden for den summeret sandsynlighed af alle alpha og beta
    % værdier ved intensitet x. 

    %p_t (lambda|x,r)
    PosteriorNextTrailSuccess = PosteriorNextTrailSuccess./DenominatorS; 
    PosteriorNextTrialFailure = PosteriorNextTrialFailure./DenominatorF;
    
    %Nødvendig i trin 4: 
    %pt(success|x) 
    pSuccessGivenx = sum(sum(sum(sum(pdf5D.*LUT,1),2),3),4); 


