
%% Subject to be stimulated setup 
clc 
clear 
close all; 

%Parameter to simulate test subject 
paramsTrue = [4, 3.05, .02, .02]; 

%%  Initial up/down 

%Purpose: To find a rough initial estimate of the threshold, which will be 
%used to define stimulation range and a prior distribution in the running 
%PSI procedure.

%Method: transformed up/down method (2 up/2 down rule)
% - Wetherill, G.B., Levitt, H., 1965. Sequential estimation of points on a
%psychometric function

%Setup 
UD.PF = @LogisticFunc;
UD.up = 2;                     %increase after 2 consecutive wrong
UD.down = 2;                   %decrease after 2 consecutive right
UD.stepSizeDown = 0.4;        
UD.stepSizeUp = 0.4;
UD.numOfStim = 20;
UD.x(1) = 1;                 %Initial intensity

while length(UD.x) <= UD.numOfStim
    UD.response(length(UD.x)) = rand(1) < UD.PF(paramsTrue, UD.x(end));   
    UD = UpdateUD(UD, UD.response(end)); %update UD structure
end

UDestimatiedThreshold = mean(UD.x(end-6:end-1));
fprintf('Up/down estimated threshold %2.1f \n',UDestimatiedThreshold)


%Create simple plot:
t = 1:length(UD.x)-1;
UDFig = figure;
plot(t,UD.x(1:end-1),'k');
hold on;
line([1 length(UD.x)-1], [UDestimatiedThreshold UDestimatiedThreshold],'linewidth', 2, 'linestyle', '--', 'color','k');

plot(t(UD.response == 1),UD.x(UD.response == 1),'ko', 'MarkerFaceColor','k');
plot(t(UD.response == 0),UD.x(UD.response == 0),'ko', 'MarkerFaceColor','w');
xlabel('Stimulation num');
ylabel('Stimulus Intensity');
uiwait(UDFig)
disp('Program execution resumed.');

%% Setup (PSI METHODE 1/3) 
NumStimulation = 300;
coef = 0.05;
windowSize = 5; 
grain     = 10; 
PM.PF = @LogisticFunc;
StimulationResolution = 4; 
PM.marginalize = [2]; 

%Stimulus values the method can select from
%PM.stimRange = (linspace(,PM.PF([paramsTrue(1) paramsTrue(2) 0 0],.99,'inverse'),StimulationResolution));
%PM.stimRange = linspace(UDestimatiedThreshold*0.5,UDestimatiedThreshold*2,StimulationResolution)
%PM.stimRange = round(linspace(PM.PF([UDestimatiedThreshold 2 0 0],.01,'inverse'),PM.PF([UDestimatiedThreshold 2 0 0],.99,'inverse'),StimulationResolution),4);
stepSizeStim = round((PM.PF([UDestimatiedThreshold 2 0 0],.99,'inverse')-PM.PF([UDestimatiedThreshold 2 0 0],.01,'inverse'))/(StimulationResolution-1),4); 
PM.stimRange = PM.PF([UDestimatiedThreshold 2 0 0],.01,'inverse'):...
    stepSizeStim:...
    PM.PF([UDestimatiedThreshold 2 0 0],.99,'inverse');

%Define parameter ranges to be included in posterior
%priorAlphaRange = linspace(PM.PF([UDestimatiedThreshold 2 0 0],.01,'inverse'),PM.PF([UDestimatiedThreshold 2 0 0],.99,'inverse'),grain);
stepSizeAlpha = (PM.PF([UDestimatiedThreshold 2 0 0],.99,'inverse',grain) -PM.PF([UDestimatiedThreshold 2 0 0],.01,'inverse'))/(grain-1);
priorAlphaRange = PM.PF([UDestimatiedThreshold 2 0 0],.01,'inverse'):...
   stepSizeAlpha: ...
   PM.PF([UDestimatiedThreshold 2 0 0],.99,'inverse',grain); 

%priorAlphaRange = linspace(PM.PF([paramsTrue(1) paramsTrue(2) 0 0],.01,'inverse'),PM.PF([paramsTrue(1) paramsTrue(2) 0 0],.99,'inverse'),grain);
priorBetaRange =  linspace(log10(.0625),log10(16),grain); %OBS. Stated in Log!
priorGammaRange = .02;  
priorLambdaRange = .02; 
[PM.priorAlphas, PM.priorBetas, PM.priorGammas, PM.priorLambdas] = ndgrid(priorAlphaRange,priorBetaRange,priorGammaRange,priorLambdaRange);

%Prior uniform 
    % "First, a prior probability distribution p0(lambda) for the 
    % psychometric functions must be set up" [Kontsevich]
    uniformPrior = ones(length(priorAlphaRange),length(priorBetaRange),length(priorGammaRange),length(priorLambdaRange));
    uniformPrior = uniformPrior./numel(uniformPrior); 
    PM.pdf = uniformPrior;  

    
%Normal distributed prior 
%     normalDistPriorAlpha = pdfNormalFunc(PM.priorAlphas, UDestimatiedThreshold, 1.88); 
%     normalDistPriorBeta = pdfNormalFunc(PM.priorBetas, log10(4.4), log10(1.88)); % Jenny values (see mail) 
%     PM.pdf = normalDistPriorAlpha .*normalDistPriorBeta; 
%     PM.pdf = PM.pdf./sum(PM.pdf, [1 2]); 

% Uncomment to see prior dist: 
% [X,Y] = meshgrid(10.^priorBetaRange, priorAlphaRange);
% contour(X,Y,PM.pdf)
% ylabel('alpha')
% xlabel('beta')
% grid 'on'
 

%LOOK UP TABEL (LUT)
    % "Second, to speed up the method, a look-up table of conditional
    % probabilities p(r|lambda,x) should be computed" [Kontsevich]
    for a = 1:length(priorAlphaRange)
        for b = 1:length(priorBetaRange) %OBS. Not calculated in log!
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
    % doPlot = input('Do not plot (0), plot threshold (1) ?: ');
    doPlot = 1 ; 
    
%% Setup plot (PSI METHODE 2/3) 

    if (doPlot) 
        figure(1) 
        xlim([1 NumStimulation])

        %set(gcf, 'Position',  [40, 3010, 1000, 600])
        xlabel('Stimulation number') 
        ylabel('Stimulus intensity')
    end 
        
                
%% Simulate data and update method (PSI METHODE 3/3) 

% first window 
[PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
[~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx, PM.marginalize);
PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
PM.x(1) = PM.xCurrent; 

while length(PM.x) <= windowSize 
    responses(length(PM.x)) = rand(1) < PM.PF([paramsTrue(1)+coef*length(PM.x) paramsTrue(2) paramsTrue(3) paramsTrue(4)], PM.xCurrent);    %simulate observer
    %  responses(length(PM.x)) = rand(1) < PM.PF([paramsTrue(1)*coef^(length(PM.x)) paramsTrue(2) paramsTrue(3) paramsTrue(4)], PM.xCurrent);    %simulate observer

    PM = UpdateFunc(PM, responses(end) );
    
    if (doPlot) 
        figure(1) 
        hold on;        
        plot(1:length(PM.x)-1, PM.threshold, 'b')
        plot(length(PM.x)-1, paramsTrue(1)+coef*length(PM.x), '.', 'color','#B1B1B1', 'linewidth',0.1) % FEJL, den er forskudt med en. 
        % plot(length(PM.x)-1, paramsTrue(1)*coef^(length(PM.x)), '.', 'color','#B1B1B1', 'linewidth',0.1) % FEJL, den er forskudt med en. 

        figure(1) 
        hold on; 
        if responses(end) 
            plot(length(PM.x)-1,PM.x(end-1),'ok','MarkerFaceColor','k');
        else 
            plot(length(PM.x)-1,PM.x(end-1),'or');
        end
    end
end 


%% Magi
%length(PM.x) = window + 1 

for curPos =1:NumStimulation-windowSize % 1:95 --> svarer til 6:100
    responses(length(PM.x)) = rand(1) < PM.PF([paramsTrue(1)+coef*length(PM.x) paramsTrue(2) paramsTrue(3) paramsTrue(4)], PM.x(end));
    % responses(length(PM.x)) = rand(1) < PM.PF([paramsTrue(1)*coef^(length(PM.x)) paramsTrue(2) paramsTrue(3) paramsTrue(4)], PM.x(end));

%Update stim range
    findOldStim = PM.stimRange(ismember(round(PM.stimRange,3),round(PM.x(curPos+1:curPos+windowSize),3) )); % 1 : 1+5-1 --> 1:5
    midStim = PM.stimRange(StimulationResolution/2); 
    if PM.threshold(end) > midStim
        while PM.threshold(end) > midStim 
            midStim = midStim+stepSizeStim;           
        end 
    else
        while PM.threshold(end) < midStim 
            midStim = midStim-stepSizeStim;           
        end
    end     
    
    newStimRange = (-StimulationResolution/2:1:StimulationResolution/2)*stepSizeStim+midStim;
    findOverlapStim = not(ismember(round(findOldStim,4), round(newStimRange,4)));
    connectedStimRange = sort([findOldStim(findOverlapStim) ,newStimRange],'ascend');
    PM.stimRange = connectedStimRange; 
    
    
    while(~isempty(find(mod(round(PM.stimRange-midStim,4), stepSizeStim)) ~= 0))
           kage = round(PM.stimRange-midStim,3); 
           for i = 1:length( kage )-1 
               stepSizes(i) = kage(1+i) - kage(i);
           end
           stepSizes
           warning('StepSize WRONG!')
           quit(code)
    end 
    
    
%Update alpha range 
    midAlpha = priorAlphaRange(grain/2); 
    stepSizeAlpha = priorAlphaRange(2) - priorAlphaRange(1); 
    
    if PM.threshold(end) > midAlpha
        while PM.threshold(end) > midAlpha 
            midAlpha = midAlpha+stepSizeAlpha;           
        end 
    else
        while PM.threshold(end) < midAlpha 
            midAlpha = midAlpha-stepSizeAlpha;          
        end
    end 
    newAlphaRange = (-grain/2:1:grain/2)*stepSizeAlpha+midAlpha;   
    findOverlapAlphaRange = not(ismember(round(priorAlphaRange,3),round(newAlphaRange,3)));  %Find value of stim need in new stim array .
    connectedAlphaRange = [priorAlphaRange(findOverlapAlphaRange), newAlphaRange];
    priorAlphaRange = sort(connectedAlphaRange,'ascend') ;

%Updade lut with new alpha 
    PM.LUT = []; 
    for a = 1:length(priorAlphaRange)
        for b = 1:length(priorBetaRange) %OBS. Not calculated in log!
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

%Reset pdf 
    uniformPrior = ones(length(priorAlphaRange),length(priorBetaRange),length(priorGammaRange),length(priorLambdaRange));
    uniformPrior = uniformPrior./numel(uniformPrior); 
    PM.pdf = uniformPrior;  
    
[PM.priorAlphas, PM.priorBetas, PM.priorGammas, PM.priorLambdas] = ndgrid(priorAlphaRange,priorBetaRange,priorGammaRange,priorLambdaRange);

    %curPos = 81 
    % PM.x(82:86) = 7.9226    7.4390    8.6480    8.4062    5.2628
    for curWinPos = curPos+1:windowSize+curPos   %  82    83    84    85    86
        [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
        
        if responses(curWinPos) == 1
            PM.pdf = PM.PosteriorNextTrailSuccess(:,:,:,:,find(round(PM.stimRange,4) == round(PM.x(curWinPos),4))); 
        else
            PM.pdf = PM.PosteriorNextTrialFailure(:,:,:,:,find(round(PM.stimRange,4) == round(PM.x(curWinPos),4)));
        end
        PM.pdf = PM.pdf./sum(sum(sum(sum(PM.pdf))));
    end 
    [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
    PM.threshold(length(PM.x)) = sum(sum(sum(sum(PM.priorAlphas.*PM.pdf))));
    
    
    [~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx, PM.marginalize);
    PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
      
    PM.x(length(PM.x)+1) = PM.xCurrent;
    
    if (doPlot) 
        figure(1) 
        hold on;        
        plot(1:length(PM.x)-1, PM.threshold, 'b')
        plot(length(PM.x)-1, paramsTrue(1)+coef*length(PM.x), '.', 'color','#B1B1B1', 'linewidth',0.1)
        % plot(length(PM.x)-1, paramsTrue(1)*coef^(length(PM.x)), '.', 'color','#B1B1B1', 'linewidth',0.1)

        figure(1) 
        hold on; 
        if responses(end) 
            plot(length(PM.x)-1,PM.x(end-1),'ok','MarkerFaceColor','k');
        else 
            plot(length(PM.x)-1,PM.x(end-1),'or');
        end
    end
end 

%% Measure residual 
% 
% s = zeros(1, length(PM.threshold)); % Preallocation
% for i = 1:length(PM.threshold)  
%     s(i) = (PM.threshold(i) - paramsTrue(1)+coef*i)^2;
%     % s1(i) = (PM.threshold(i) - paramsTrue(1)*(coef)^i)^2;
% end 
% 
% fprintf('Residual %2.2f. \n',sum(s))

%% Maximum likelihood (1/2 - bruteforce search) 
stimIntst = PM.x(1:end-1); 
%NumPos = [59 53 68 83 92 99]; %number of trials with correct response.
%OutOfNum = [100 100 100 100 100 100]; %number of stimulus in each stimLevel

%Define sg for brute force search 
sg.alpha = linspace(4,7,100);
sg.beta = 3.05; %linspace(log10(.0625),log10(5),100); 
sg.delta = linspace(0,1,100); %habituation factor 
sg.gamma = 0.02;
sg.lambda = 0.02;

% Find likelihood grid 
likelihoodGrid = zeros(length(sg.alpha),length(sg.delta)); 

for a = 1:length(sg.alpha)
    for d = 1:length(sg.delta)
        temp = zeros(1,length(stimIntst)); % Preallocation
        for stimNum = 1:length(stimIntst)
            if responses(stimNum) == true
                temp(stimNum) = log(PM.PF([sg.alpha(a)+sg.delta(d)*stimNum, sg.beta, sg.gamma, sg.lambda], stimIntst(stimNum)));
            else
                temp(stimNum) = log(1-PM.PF([sg.alpha(a)+sg.delta(d)*stimNum, sg.beta, sg.gamma, sg.lambda], stimIntst(stimNum)));
            end 
        end
        likelihoodGrid(a,d)= sum(temp); 
    end
end 
[M,I] = max(likelihoodGrid(:)); 
[I_row, I_col] = ind2sub(size(likelihoodGrid),I); 
parameters = [sg.alpha(I_row) sg.delta(I_col)]; 

%estimatedThreshold = parameters(1)+(1:length(PM.x)-1)*parameters(2);

figure(1)
hold on; 
plot(1:(length(PM.x)-1), parameters(1)+(1:length(PM.x)-1)*parameters(2), '.', 'color','r', 'linewidth',0.1) % FEJL, den er forskudt med en. 



