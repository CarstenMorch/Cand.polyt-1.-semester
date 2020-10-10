
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

%% Look up tabel  p(success if lampda given x)
clc 
clear 
disp(' >> Look up tabel << '); 

grain = 201; 
PF = @PAL_Gumbel;

%Stimulus values the method can select from
PM.stimRange = (linspace(PF([0 1 0 0],.1,'inverse'),PF([0 1 0 0],.9999,'inverse'),21));

%Define parameter ranges to be included in posterior
PM.priorAlphaRange = linspace(PF([0 1 0 0],.1,'inverse'),PF([0 1 0 0],.9999,'inverse'),grain);
PM.priorBetaRange =  linspace(log10(.0625),log10(16),grain); %OBS ANGIVET I LOG!
PM.priorGammaRange = 0.5;  
PM.priorLambdaRange = .02; 
PM.gammaEQlambda = 0; 

%parameter to simulate observer
paramsGen = [0, 1, .5, .02]; 

%PDF :
    PM.prior = ones(length(PM.priorAlphaRange),length(PM.priorBetaRange),length(PM.priorGammaRange),length(PM.priorLambdaRange));
    PM.prior = PM.prior./numel(PM.prior); % number of elements
    PM.pdf = PM.prior; 

 
%LOOK UP TABEL   : 
    for a = 1:length(PM.priorAlphaRange)
        for b = 1:length(PM.priorBetaRange) %OBS. Udregnes ikke i log!
            for g = 1:length(PM.priorGammaRange)
                for L = 1:length(PM.priorLambdaRange) 
                    for sLevel = 1:length(PM.stimRange)
                        %sandsyndligheden for at for korrekt response ved en given parameter sammensætning og intensitet.  
                        LUT(a,b,g,L,sLevel) = PF([PM.priorAlphaRange(a), 10.^PM.priorBetaRange(b), PM.priorGammaRange(g), PM.priorLambdaRange(L)], PM.stimRange(sLevel));
                    end
                end
            end
        end 
    end     
    
%% Trin 1 og 2 (posterior probability)
    pdf5D = repmat(PM.pdf, [1 1 1 1 length(PM.stimRange)]); % 201 201 --> 201 201 1 1 21 
    PosteriorTplus1givenSuccess = pdf5D.*LUT; % p(succes/lambda, x)
    PosteriorTplus1givenFailure = pdf5D-PosteriorTplus1givenSuccess;
   
    size(PosteriorTplus1givenSuccess)
    
    Denominator = squeeze(sum(sum(sum(sum(PosteriorTplus1givenSuccess,1),2),3),4)); %sandsynligheden for x ift all parameter, return the 1,1,1,1,21 --> 21,1
    Denominator = repmat(Denominator, [1 size(pdf5D,1) size(pdf5D,2) size(pdf5D,3) size(pdf5D,4)]); %[21,1] --> [21, 201, 201]    
    Denominator = permute(Denominator, [2 3 4 5 1]); %[21, 201, 201] --> [201 201 1 1 21]
    
% Sandsynligheden for alpha og beta ved intensitet x divideret med 
% sandsynligheden for den summeret sandsynlighed af alle alpha og beta
% værdier ved intensitet x. 
% p_t (lambda|x,r)
    PosteriorTplus1givenSuccess = PosteriorTplus1givenSuccess./Denominator; 
    
    Denominator = squeeze(sum(sum(sum(sum(PosteriorTplus1givenFailure,1),2),3),4));
    Denominator = repmat(Denominator, [1 size(pdf5D,1) size(pdf5D,2) size(pdf5D,3) size(pdf5D,4)]);
    Denominator = permute(Denominator, [2 3 4 5 1]);
    PosteriorTplus1givenFailure = PosteriorTplus1givenFailure./Denominator;
    
    %pt(success|x)
    pSuccessGivenx = sum(sum(sum(sum(pdf5D.*LUT,1),2),3),4); 
    

%% Trin 3,4 and 5 (ENTROPY)
    %success
    EntropyS = PosteriorTplus1givenSuccess.*log(PosteriorTplus1givenSuccess);
    EntropyS(isnan(EntropyS)) = 0;          %effectively defines 0.*log(0) to equal 0.
    EntropyS = sum(EntropyS,'all');
    EntropyS = -EntropyS;
    
    %Failure
    EntropyF = PosteriorTplus1givenFailure.*log(PosteriorTplus1givenFailure);
    EntropyF(isnan(EntropyF)) = 0;          %effectively defines 0.*log(0) to equal 0.
    EntropyF = sum(EntropyF,'all');
    EntropyF = -EntropyF;
    
% E[Ht(x)]: 
    Entropy = EntropyS.*pSuccessGivenx + EntropyF.*(1-pSuccessGivenx); % [1,1,1,1,21]
    
% Find minimum: 
    [minEntropy, newIntensityLevel] = min(squeeze(Entropy)); % minentropy: the value, i; index position. 
    PM.xCurrent = PM.stimRange(newIntensityLevel);
    PM.x(1) = PM.xCurrent; %The first stimulation value!! what up!

    response = rand(1) < PF(paramsGen, PM.xCurrent);    %simulate observer

