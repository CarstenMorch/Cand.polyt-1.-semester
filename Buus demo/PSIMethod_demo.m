%% Setup (PSI METHODE 1/3) 
clc 
clear 
close all; 

NumTrials = 30;
NumTrials2 = 60;
NumTrials3 = 90;
NumTrials4 = 120;
NumTrials5 = 150;
NumTrials6 = 180;

habituationsConstant = 0.02;

grain     = 50; 
PM.PF = @LogisticFunc;
StimulationResolution = 50; 

%parameter to simulate observer
paramsGen = [10, 2, .02, .02]; 

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
    clear grain StimulationResolution 
    
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
        plot([NumTrials NumTrials], [min(PM.stimRange) max(PM.stimRange)],'r')
        hold on
        plot([NumTrials2 NumTrials2], [min(PM.stimRange) max(PM.stimRange)],'r')       
        hold on
        plot([NumTrials3 NumTrials3], [min(PM.stimRange) max(PM.stimRange)],'r')
        xlim([1 NumTrials6])
        ylim([min(PM.stimRange) max(PM.stimRange)])
        xlabel('Trial number') 
        ylabel('threshold estimation')

        figure(2) 
        [X,Y] = meshgrid(priorAlphaRange, priorBetaRange);
     
        
        figure(3) 
        [X,Y] = meshgrid(priorAlphaRange, 10.^priorBetaRange);
        StimLevelsFine = [min(PM.stimRange):(max(PM.stimRange) - min(PM.stimRange))./1000:max(PM.stimRange)];
        Fit_correct = PM.PF(paramsGen, StimLevelsFine); 
    end 
                  
%% Simulate data and update method (PSI METHODE 3/3) 


while length(PM.x) <= NumTrials
    % Simulate habituation with 0.05
    response = rand(1) < PM.PF([10 + habituationsConstant.*length(PM.x), 2, .02, .02], PM.xCurrent);    %simulate observer
    
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

disp('finish first trail')
disp('Threshold estimate first trial:')
threshold1 = PM.threshold(end)

disp('Slope estimate :')
10.^PM.slope(end)           %PM.slope is in log10 units of beta parameter

%% Estimates threshold and slope for second trial
clear PM.x
PM.pdf = prior;
[PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
% Step 3 to 5 
[~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
% next stimulation intensity
PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
PM.x(1) = PM.xCurrent;


for i=NumTrials:NumTrials2
 
    % Simulate habituation with 0.05
    response = rand(1) < PM.PF([10 + habituationsConstant.*(length(PM.x)), 2, .02, .02], PM.xCurrent);    %simulate observer
    
    %update PM based on response
    PM = UpdateFunc(PM, response); 
    % Test --> PMtest = PAL_AMPM_updatePM(PMtest,response);
    responses(length(PM.x)-1) = response; 
    
    % plot? 
    if (doPlot) 
        figure(1) 
        hold on; 
        plot([NumTrials+1:length(PM.x)-1],PM.threshold(NumTrials+1:end), 'b')
          
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

disp('finish trial 2 ')
disp('Threshold estimate:')
threshold2 = PM.threshold(end)

disp('Slope estimate trial 2 :')
10.^PM.slope(end)           %PM.slope is in log10 units of beta parameter

%% Estimate threshold and slope for third trial
PM.pdf = prior;
[PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
% Step 3 to 5 
[~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
% next stimulation intensity
PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
PM.x(1) = PM.xCurrent;


for i=NumTrials2:NumTrials3
 
    % Simulate habituation with 0.05
    response = rand(1) < PM.PF([10 + habituationsConstant.*(length(PM.x)), 2, .02, .02], PM.xCurrent);    %simulate observer
    
    %update PM based on response
    PM = UpdateFunc(PM, response); 
    % Test --> PMtest = PAL_AMPM_updatePM(PMtest,response);
    responses(length(PM.x)-1) = response; 
    
    % plot? 
    if (doPlot) 
        figure(1) 
        hold on; 
        plot([NumTrials2+1:length(PM.x)-1],PM.threshold(NumTrials2+1:end), 'b')
        %plot(PM.threshold, 'b')
          
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

disp('finish trial 3 ')
disp('Threshold estimate:')
threshold3 = PM.threshold(end)

disp('Slope estimate trial 3:')
10.^PM.slope(end)           %PM.slope is in log10 units of beta parameter

%% Estimate threshold and slope for fourth trial
PM.pdf = prior;
[PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
% Step 3 to 5 
[~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
% next stimulation intensity
PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
PM.x(1) = PM.xCurrent;


for i=NumTrials3:NumTrials4
 
    % Simulate habituation with 0.05
    response = rand(1) < PM.PF([10 + habituationsConstant.*(length(PM.x)), 2, .02, .02], PM.xCurrent);    %simulate observer
    
    %update PM based on response
    PM = UpdateFunc(PM, response); 
    % Test --> PMtest = PAL_AMPM_updatePM(PMtest,response);
    responses(length(PM.x)-1) = response; 
    
    % plot? 
    if (doPlot) 
        figure(1) 
        hold on; 
        plot([NumTrials3+1:length(PM.x)-1],PM.threshold(NumTrials3+1:end), 'b')
        %plot(PM.threshold, 'b')
          
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

disp('finish trial 4 ')
disp('Threshold estimate:')
threshold4 = PM.threshold(end)

disp('Slope estimate trial 4:')
10.^PM.slope(end)           %PM.slope is in log10 units of beta parameter

%% Estimate threshold and slope for fifth trial
PM.pdf = prior;
[PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
% Step 3 to 5 
[~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
% next stimulation intensity
PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
PM.x(1) = PM.xCurrent;


for i=NumTrials4:NumTrials5
 
    % Simulate habituation with 0.05
    response = rand(1) < PM.PF([10 + habituationsConstant.*(length(PM.x)), 2, .02, .02], PM.xCurrent);    %simulate observer
    
    %update PM based on response
    PM = UpdateFunc(PM, response); 
    % Test --> PMtest = PAL_AMPM_updatePM(PMtest,response);
    responses(length(PM.x)-1) = response; 
    
    % plot? 
    if (doPlot) 
        figure(1) 
        hold on; 
        plot([NumTrials4+1:length(PM.x)-1],PM.threshold(NumTrials4+1:end), 'b')
        %plot(PM.threshold, 'b')
          
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

disp('finish trial 5 ')
disp('Threshold estimate:')
threshold5 = PM.threshold(end)

disp('Slope estimate trial 5:')
10.^PM.slope(end)           %PM.slope is in log10 units of beta parameter

%% Estimate threshold and slope for six trial
PM.pdf = prior;
[PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
% Step 3 to 5 
[~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
% next stimulation intensity
PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
PM.x(1) = PM.xCurrent;


for i=NumTrials5:NumTrials6
 
    % Simulate habituation with 0.05
    response = rand(1) < PM.PF([10 + habituationsConstant.*(length(PM.x)), 2, .02, .02], PM.xCurrent);    %simulate observer
    
    %update PM based on response
    PM = UpdateFunc(PM, response); 
    % Test --> PMtest = PAL_AMPM_updatePM(PMtest,response);
    responses(length(PM.x)-1) = response; 
    
    % plot? 
    if (doPlot) 
        figure(1) 
        hold on; 
        plot([NumTrials5+1:length(PM.x)-1],PM.threshold(NumTrials5+1:end), 'b')
        %plot(PM.threshold, 'b')
          
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

disp('finish trial 6 ')
disp('Threshold estimate:')
threshold6 = PM.threshold(end)

disp('Slope estimate trial 6:')
10.^PM.slope(end)           %PM.slope is in log10 units of beta parameter

threshold = [threshold1;threshold2;threshold3;threshold4;threshold5;threshold6]

plot(threshold,'o')
xlabel('Trial number');
ylabel('Threshold value');
title 'Threshold estimation for habituation'
% Test :
% disp('Threshold estimate - palamedes:')
% PMtest.threshold(end)
% 
% disp('Slope estimate - palamedes:')
% 10.^PMtest.slope(end)           %PM.slope is in log10 units of beta parameter
