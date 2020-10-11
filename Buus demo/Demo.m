
%% Create a Logistic Function 
clc;
clear;

StimLevels = [0:0.5:10]; 
coef = [3 1 0.1 0]; 
pcorrect = LogisticFunc(coef, StimLevels); 

figure(1)
plot(StimLevels, pcorrect, 'ko')

%% Brute force search (Maximum likelihood 1/2)
clc 
clear 

PF = @LogisticFunc;

StimLevels = [.01 .03 .05 .07 .09 .11];
NumPos = [59 53 68 83 92 99]; %number of trials with correct response.
OutOfNum = [100 100 100 100 100 100]; %number of stimulus in each stimLevel

%Define searchgrid for brute force search 
searchGrid.alpha = [0.01:0.001:0.11];
searchGrid.beta = logspace(0,3,101);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0.02;

[paramsValues] = BruteForceSearchFunc(StimLevels, NumPos, OutOfNum, ...
    searchGrid, PF)
  
%% Steepest ascent (Maximum likelihood 2/2)

SteepestAscentFunc(paramsValues, StimLevels, NumPos, OutOfNum, PF)

%% PLot data. 

PropCorrectData = NumPos./OutOfNum;
StimLevelsFine = [min(StimLevels):(max(StimLevels) - min(StimLevels))./1000:max(StimLevels)];
Fit = PF(paramsValues, StimLevelsFine);
plot(StimLevelsFine,Fit,'g-','linewidth',2);
hold on;
plot(StimLevels, PropCorrectData,'k.','markersize',40);

%% Buffer zone 




%% Setup
clc 
clear 

grain = 201; 
PM.PF = @PAL_Gumbel;

%Stimulus values the method can select from
PM.stimRange = (linspace(PM.PF([0 1 0 0],.1,'inverse'),PM.PF([0 1 0 0],.9999,'inverse'),21));

%Define parameter ranges to be included in posterior
priorAlphaRange = linspace(PM.PF([0 1 0 0],.1,'inverse'),PM.PF([0 1 0 0],.9999,'inverse'),grain);
priorBetaRange =  linspace(log10(.0625),log10(16),grain); %OBS ANGIVET I LOG!
priorGammaRange = 0.5;  
priorLambdaRange = .02; 
gammaEQlambda = 0; 

%parameter to simulate observer
paramsGen = [0, 1, .5, .02]; 

%PDF
    %First, a prior probability distribution p0(lambda) for the psychometric 
    %functions must be set up.
    prior = ones(length(priorAlphaRange),length(priorBetaRange),length(priorGammaRange),length(priorLambdaRange));
    prior = prior./numel(prior); % number of elements
    PM.pdf = prior;  
    
    
    
%LOOK UP TABEL 
    %Second, to speed up the method, a look-up table of conditional
    %probabilities p(r|lambda,x) should be computed
    for a = 1:length(priorAlphaRange)
        for b = 1:length(priorBetaRange) %OBS. Udregnes ikke i log!
            for g = 1:length(priorGammaRange)
                for L = 1:length(priorLambdaRange) 
                    for sLevel = 1:length(PM.stimRange)
                        %sandsyndligheden for at for korrekt response ved en given parameter sammensætning og intensitet.  
                        PM.LUT(a,b,g,L,sLevel) = PM.PF([priorAlphaRange(a), 10.^priorBetaRange(b), priorGammaRange(g), priorLambdaRange(L)], PM.stimRange(sLevel));
                    end
                end
            end
        end 
    end
    clear a b g L sLevel 
    
%% Fist time 
    % trin 1 til 2 
    [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
    % trin 3 til 5 
    [minEntropy, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
    % næste stimulering intensitet
    PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
    PM.x(1) = PM.xCurrent;

    clear newIntensityIndexPosition minEntropy
    
%% Simulate data 
NumTrials = 20;

while length(PM.x) < NumTrials
    response = rand(1) < PM.PF(paramsGen, PM.xCurrent);    %simulate observer

    %update PM based on response
    PM = noNameFunc(PM, response); 
    %PM = PAL_AMPM_updatePM(PM,response);
    break; 
end

    

