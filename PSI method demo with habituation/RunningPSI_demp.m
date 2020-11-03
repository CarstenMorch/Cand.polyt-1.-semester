
%% Setup (PSI METHODE 1/3) 
clc 
clear 
close all; 

disp('Demo shows PSI method')
disp('Simulated data habituate between stimulation with 0.1% factor')
disp('Habituation continuous between trials'); 

NumStimulation = 50;
coef = 1.002;
windowSize = 20; 
grain     = 50; 
PM.PF = @LogisticFunc;
StimulationResolution = 50; 

%parameter to simulate observer
paramsGen = [5, 2, .02, .02]; 

%Stimulus values the method can select from
PM.stimRange = (linspace(PM.PF([paramsGen(1) paramsGen(2) 0 0],.01,'inverse'),PM.PF([paramsGen(1) paramsGen(2) 0 0],.99,'inverse'),StimulationResolution));

%Define parameter ranges to be included in posterior
priorAlphaRange = linspace(PM.PF([paramsGen(1) paramsGen(2) 0 0],.01,'inverse'),PM.PF([paramsGen(1) paramsGen(2) 0 0],.99,'inverse'),grain);
priorBetaRange =  linspace(log10(.0625),log10(5),grain); %OBS. Stated in Log!
priorGammaRange = .02;  
priorLambdaRange = .02; 

[PM.priorAlphas, PM.priorBetas, PM.priorGammas, PM.priorLambdas] = ndgrid(priorAlphaRange,priorBetaRange,priorGammaRange,priorLambdaRange);

%PDF
    % "First, a prior probability distribution p0(lambda) for the 
    % psychometric functions must be set up" [Kontsevich]
    prior = ones(length(priorAlphaRange),length(priorBetaRange),length(priorGammaRange),length(priorLambdaRange));
    prior = prior./numel(prior); 
    PM.pdf = prior;  
  
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
    clear StimulationResolution  
    doPlot = input('Do not plot (0), plot threshold (1) ?: ');

    
%% Setup plot (PSI METHODE 2/3) 

    if (doPlot) 
        figure(1) 
        xlim([1 NumStimulation])
        ylim([min(PM.stimRange) max(PM.stimRange)])

        %set(gcf, 'Position',  [40, 3010, 1000, 600])
        xlabel('Stimulation number') 
        ylabel('Stimulus intensity')
    end 
                  
                
%% Simulate data and update method (PSI METHODE 3/3) 

% first window 
[PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
[~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
PM.x(1) = PM.xCurrent; 

while length(PM.x) <= windowSize 
    responses(length(PM.x)) = rand(1) < PM.PF([paramsGen(1)*coef^(length(PM.x)) paramsGen(2) paramsGen(3) paramsGen(4)], PM.xCurrent);    %simulate observer
    PM = UpdateFunc(PM, responses(end) );
    
    if (doPlot) 
        figure(1) 
        hold on;        
        plot(1:length(PM.x)-1, PM.threshold, 'b')
        plot(length(PM.x)-1, paramsGen(1)*coef^(length(PM.x)), '.', 'color','#B1B1B1', 'linewidth',0.1)

        figure(1) 
        hold on; 
        if responses(end) 
            plot(length(PM.x)-1,PM.x(end-1),'ok','MarkerFaceColor','k');
        else 
            plot(length(PM.x)-1,PM.x(end-1),'or');
        end
    end
end 


for curPos =1:NumStimulation-windowSize 
    responses(length(PM.x)) = rand(1) < PM.PF([paramsGen(1)*coef^(length(PM.x)) paramsGen(2) paramsGen(3) paramsGen(4)], PM.x(end));    

    PM.pdf = prior;  
    
    for curWinPos = curPos+1:windowSize+curPos   
        [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
        if responses(curWinPos) == 1
            PM.pdf = PM.PosteriorNextTrailSuccess(:,:,:,:,find(PM.stimRange == PM.x(curWinPos))); 
        else
            PM.pdf = PM.PosteriorNextTrialFailure(:,:,:,:,find(PM.stimRange == PM.x(curWinPos)));
        end
        PM.pdf = PM.pdf./sum(sum(sum(sum(PM.pdf)))); 
    end 
    
    [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
    PM.threshold(length(PM.x)) = sum(sum(sum(sum(PM.priorAlphas.*PM.pdf))));

    [~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
    PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
      
    PM.x(length(PM.x)+1) = PM.xCurrent;
    
    if (doPlot) 
        figure(1) 
        hold on;        
        plot(1:length(PM.x)-1, PM.threshold, 'b')
        plot(length(PM.x)-1, paramsGen(1)*coef^(length(PM.x)), '.', 'color','#B1B1B1', 'linewidth',0.1)

        figure(1) 
        hold on; 
        if responses(end) 
            plot(length(PM.x)-1,PM.x(end-1),'ok','MarkerFaceColor','k');
        else 
            plot(length(PM.x)-1,PM.x(end-1),'or');
        end
    end
end 

