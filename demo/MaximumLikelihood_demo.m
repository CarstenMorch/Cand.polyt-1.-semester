
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