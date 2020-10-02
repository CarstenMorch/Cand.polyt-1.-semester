

%% Create a Logistic Function 
clc;
clear;

StimLevels = [0:0.5:10]; 
pcorrect = LogisticFunction([3 1 0.1 0], StimLevels); 

figure(1)
plot(StimLevels, pcorrect, 'ko')

%% Maximum likelihood (Palamedes) 
clc 
clear

StimLevels = [.01 .03 .05 .07 .09 .11];
NumPos = [59 53 68 83 92 99]; %number of trials with correct response.
OutOfNum = [100 100 100 100 100 100]; %number of stimulus in each stimLevel

PF = @PAL_Logistic;
PD = @LogisticFunction;

paramsFree = [1 1 0 0]; %determine free or fixed parameters

searchGrid.alpha = [0.01:0.001:0.11];
searchGrid.beta = logspace(0,3,101);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0.02;

[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum,searchGrid, paramsFree,PF)

PropCorrectData = NumPos./OutOfNum;
StimLevelsFine = [min(StimLevels):(max(StimLevels) - min(StimLevels))./1000:max(StimLevels)];
Fit = PF(paramsValues, StimLevelsFine);
plot(StimLevelsFine,Fit,'g-','linewidth',2);
hold on;
plot(StimLevels, PropCorrectData,'k.','markersize',40);
set(gca, 'fontsize',12);
axis([0 .12 .4 1]);

%% Maximum likelihood function (homemade) 
%Se article 'A simplex method for function minimization' 

clc 
clear 

    PF = @PAL_Logistic;
%     StimLevels = [-3:1:3];
%     NumPos = [55 55 66 75 91 94 97];    %observer data
%     OutOfNum = 100.*ones(size(StimLevels));
    
    StimLevels = [.01 .03 .05 .07 .09 .11];
    NumPos = [59 53 68 83 92 99]; %number of trials with correct response.
    OutOfNum = [100 100 100 100 100 100]; %number of stimulus in each stimLevel

    searchGrid.alpha = [0.01:0.001:0.11];
    searchGrid.beta = logspace(0,3,101);
    searchGrid.gamma = 0.5;
    searchGrid.lambda = 0.02;
   
    
   disp('palamedes') 
   [paramsValues] = PAL_PFML_BruteForceFit(StimLevels, NumPos, ...
   OutOfNum, searchGrid, PF)

   disp('homemade') 
   [paramsValues] = BruteForceSearchFunction(StimLevels, NumPos, ...
       OutOfNum, searchGrid, PF)
  
%% steepest ascent 

% ALGORITHM PARAMETERS 
dx    = 0.001; 
dy    = 0.001; 
dz    = 0.001;
dm    = 0.001; 
alpha = 0.1; 

% Initial Guess 
x0 = paramsValues(1);
y0 = paramsValues(2);
z0 = paramsValues(3); 
m0 = paramsValues(4);


% PREFROM ALGORITHM
tol = 1e-100; 
g   = [inf;inf;inf];

n = 0;
while norm(g) > tol 
    % compute gradient 
    f1 = SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0-dx/2, y0, z0, m0], PF); 
    f2 = SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0+dx/2, y0, z0, m0], PF);
    gx = (f2 - f1)/dx; 
    
    f1 = SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0, y0-dy/2, z0, m0], PF); 
    f2 = SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0, y0+dy/2, z0, m0], PF);
    gy = (f2 - f1)/dy; 
    
    f1 = SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0,y0, z0-dz/2, m0], PF); 
    f2 = SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0,y0, z0+dz/2, m0], PF);
    gz = (f2 - f1)/dz; 
    
    f1 = SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0,y0, z0, m0-dm/2], PF); 
    f2 = SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0,y0, z0, m0-dm/2], PF);
    gm = (f2 - f1)/dm;
    
    g = [ gx ; gy ; gz ; gm];
   
    
    % update Postion of Guess
    x0 = x0 + alpha*gx; 
    y0 = y0 + alpha*gy;     
    z0 = z0 + alpha*gz; 
    m0 = m0 + alpha*gm; 

    % counter
    n = n+1; 
end 

% REPORT THAN ANSWER
[x0 y0 z0 m0]
SingleValueLikelihoodFunction(StimLevels,NumPos, OutOfNum,[x0,y0, z0, m0], PF)
C = {'number of iteratation :', n };
sprintf('%s %d ',C{:})  



%%

PropCorrectData = NumPos./OutOfNum;
StimLevelsFine = [min(StimLevels):(max(StimLevels) - min(StimLevels))./1000:max(StimLevels)];
Fit = PF(paramsValues, StimLevelsFine);
plot(StimLevelsFine,Fit,'g-','linewidth',2);
hold on;
plot(StimLevels, PropCorrectData,'k.','markersize',40);
set(gca, 'fontsize',12);
axis([0 .12 .4 1]);






