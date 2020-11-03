
%% Setup (PSI METHODE 1/3) 
clc 
clear 
close all; 

NumTrials = 100;
grain     = 50; 
PM.PF = @LogisticFunc;
StimulationResolution = 50; 

%parameter to simulate observer
paramsGen = [1, 2, .05, .02];  %[alpha, beta, gamma, lambda]

%Stimulus values the method can select from
PM.stimRange = (linspace(PM.PF([paramsGen(1) paramsGen(2) 0 0],.1,'inverse'),PM.PF([paramsGen(1) paramsGen(2) 0 0],.9999,'inverse'),StimulationResolution));

%Define parameter ranges to be included in posterior
priorAlphaRange = linspace(PM.PF([paramsGen(1) paramsGen(2) 0 0],.1,'inverse'),PM.PF([paramsGen(1) paramsGen(2) 0 0],.9999,'inverse'),grain);
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
    

     %TEST: 
%     PMtest = PAL_AMPM_setupPM('priorAlphaRange',priorAlphaRange,...
%                   'priorBetaRange',priorBetaRange,...
%                   'priorGammaRange',priorGammaRange,...
%                   'priorLambdaRange',priorLambdaRange,...
%                   'numtrials',NumTrials,...
%                   'PF' , PM.PF,...
%                   'stimRange',PM.stimRange);  
    
    clear a b g L sLevel 
    clear grain StimulationResolution prior 
    
    doPlot = input('Do not plot (0), plot threshold (1)?: ');

    
%% First time (PSI METHODE 2/3) 
    % Step 1 to 2 
    [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
    % Step 3 to 5 
    [~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
    % next stimulation intensity
    PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
    PM.x(1) = PM.xCurrent;

    clear newIntensityIndexPosition 
    
    if (doPlot) 
        figure(1) 
        xlim([1 NumTrials])
        ylim([min(PM.stimRange) max(PM.stimRange)])
        set(gcf, 'Position',  [40, 300, 1000, 600])
        xlabel('Trial number') 
        ylabel('threshold estimation')

        figure(2) 
        [X,Y] = meshgrid(priorAlphaRange, priorBetaRange);
        set(gcf, 'Position',  [1100, 580, 800, 400])
        
        figure(3) 
        [X,Y] = meshgrid(priorAlphaRange, 10.^priorBetaRange);
        StimLevelsFine = [min(PM.stimRange):(max(PM.stimRange) - min(PM.stimRange))./1000:max(PM.stimRange)];
        Fit_correct = PM.PF(paramsGen, StimLevelsFine); 
        set(gcf, 'Position',  [1100, 50, 800, 400])
    end 
                  
%% Simulate data and update method (PSI METHODE 3/3) 

while length(PM.x) <= NumTrials
    response = rand(1) < PM.PF(paramsGen, PM.xCurrent);    %simulate observer
    
    %update PM based on response
    PM = UpdateFunc(PM, response); 
    % Test --> PMtest = PAL_AMPM_updatePM(PMtest,response);
    responses(length(PM.x)-1) = response; 
    
    % plot? 
    if (doPlot) 
        figure(1) 
        hold on; 
        plot(PM.threshold, 'b')
          
        figure(1) 
        hold on; 
        if response 
            plot(length(PM.x),PM.x(end),'ok','MarkerFaceColor','k');
        else 
            plot(length(PM.x),PM.x(end),'ok');
        end
        
        figure(2) 
        contour(X,Y, flip(flip(PM.pdf),2))        
        xlabel('threshold')
        ylabel('slope')
        title('Probability density function (PDF)')
        grid on; 
               
        figure(3) 
        Fit = PM.PF([PM.threshold(end), 10.^PM.slope(end), 0.02, 0.02], StimLevelsFine);
        clf
        hold on; 
        plot(StimLevelsFine,Fit_correct,'k-','linewidth',2);
        plot(StimLevelsFine,Fit,'g-','linewidth',2);
        title('Psychophysical Function')

        for i = 2:length(responses)
            if responses(i) 
                dotSize = sum(PM.x(i) == PM.x(i:-1:1) & responses(i) == responses(i:-1:1));                                          
                plot(PM.x(i),responses(i),'ok','MarkerFaceColor','k', 'linewidth',dotSize);
            else
                dotSize = sum(PM.x(i) == PM.x(i:-1:1) & responses(i) == responses(i:-1:1));                                          
                plot(PM.x(i),responses(i),'ok','linewidth',dotSize);
            end 
        end 
        drawnow
    end
end
 
disp('finish')
disp('Threshold estimate - homemade:')
PM.threshold(end)

disp('Slope estimate - homemade :')
10.^PM.slope(end)           %PM.slope is in log10 units of beta parameter

% Test :
% disp('Threshold estimate - palamedes:')
% PMtest.threshold(end)
% 
% disp('Slope estimate - palamedes:')
% 10.^PMtest.slope(end)           %PM.slope is in log10 units of beta parameter

