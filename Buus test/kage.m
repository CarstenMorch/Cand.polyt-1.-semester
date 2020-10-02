

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
    StimLevels = [-3:1:3];
    NumPos = [55 55 66 75 91 94 97];    %observer data
    OutOfNum = 100.*ones(size(StimLevels));

    searchGrid.alpha = [-1:.1:1];
    searchGrid.beta = 10.^[-1:.1:2];
    searchGrid.gamma = .5;
    searchGrid.lambda = [0:.005:.06];
   
    
   disp('palamedes') 
   [paramsValues] = PAL_PFML_BruteForceFit(StimLevels, NumPos, ...
   OutOfNum, searchGrid, PF)

   disp('homemade') 
   [paramsValues] = BruteForceSearchFunction(StimLevels, NumPos, ...
       OutOfNum, searchGrid, PF)
   
   %%
   
%     for a = 1:length(searchGrid.alpha)
%         for b = 1:length(searchGrid.beta)
%             for g = 1:length(searchGrid.gamma)
%                 for L = 1:length(searchGrid.lambda) 
%                     parameterGrid(b,a,g,L) = {[searchGrid.alpha(a) searchGrid.beta(b) searchGrid.gamma(g)  searchGrid.lambda(L)]};
%                 end
%             end
%         end 
%     end 
   
    likelihoodGrid = zeros(length(searchGrid.alpha),length(searchGrid.beta),length(searchGrid.gamma),length(searchGrid.lambda));
    
    for a = 1:length(searchGrid.alpha)
        for b = 1:length(searchGrid.beta)
            for g = 1:length(searchGrid.gamma)
                for L = 1:length(searchGrid.lambda) 
                    for sLevel = 1:length(StimLevels)
                        temp(sLevel) = (PF([searchGrid.alpha(a), searchGrid.beta(b), searchGrid.gamma(g), searchGrid.lambda(L)], StimLevels(sLevel))^NumPos(sLevel))*((1-PF([searchGrid.alpha(a), searchGrid.beta(b), searchGrid.gamma(g), searchGrid.lambda(L)], StimLevels(sLevel)))^(OutOfNum(sLevel)-NumPos(sLevel)));
                    end
                    likelihoodGrid(a,b,g,L)=prod(temp); 
                    clear temp; 
                end
            end
        end 
    end 
    
    Amax = max(likelihoodGrid(:));                                           % Maximum Value
    Idx = find(likelihoodGrid(:) == Amax);                                   % Returns Linear Indices To All ‘Amax’ Values
    [r,c,dim3,dim4] = ind2sub(size(likelihoodGrid), Idx); 
    
    parameter = [searchGrid.alpha(r) searchGrid.beta(c) searchGrid.gamma(dim3) searchGrid.lambda(dim4)]
    
    
%   [paramsValues LL] = PAL_PFML_BruteForceFit(StimLevels, NumPos, ...
%       OutOfNum, searchGrid, PF)
  
  
  
  
  
  
  
  %%
  
  [parameterGrid.alpha, parameterGrid.beta, parameterGrid.gamma, parameterGrid.lambda] = ndgrid(searchGrid.alpha,searchGrid.beta,searchGrid.gamma,searchGrid.lambda)
  [paramsGrid.alpha, paramsGrid.beta, paramsGrid.gamma, paramsGrid.lambda] = ndgrid(searchGrid.alpha,searchGrid.beta,searchGrid.gamma,searchGrid.lambda)

  % [x, y, z, d] = ndgrid: X columns (x,:) , Y rovs (:,y), Z (:,:,z,:), D (:,:,:,d)
 

   
      
  matrix(1,1,:) = vector1;
  
    
  
  %Fit data:
  %PF(paramsGrid,StimLevels(1))
  LLspace = zeros(size(paramsGrid.alpha,1),size(paramsGrid.alpha,2),size(paramsGrid.alpha,3),size(paramsGrid.alpha,4));

              
  for level = 1:length(StimLevels)
        LLspace = LLspace + NumPos(level).*log(PF(paramsGrid,StimLevels(level)))+(OutOfNum(level)-NumPos(level)).*log(1-PF(paramsGrid,StimLevels(level)));
    end
    
    [maxim, I] = PAL_findMax(LLspace)


  [paramsValues LL] = PAL_PFML_BruteForceFit(StimLevels, NumPos, ...
      OutOfNum, searchGrid, PF)
  
  
  


%%
clc 
clear 

x = zeros(5,5,5)
x(1,2,3) = 12
[maxim I] = PAL_findMax(x)

LLspace = LLspace + NumPos(level).*log(PF(paramsGrid,StimLevels(level)))+(OutOfNum(level)-NumPos(level)).*log(1-PF(paramsGrid,StimLevels(level)));
% [x, y, z, d] = ndgrid([1 3],[10 2 2], [5 1], [20 40 30 20])
 
parameterGrid











