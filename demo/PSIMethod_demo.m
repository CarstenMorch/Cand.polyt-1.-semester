% 
% %% Setup (PSI METHODE 1/3) 
% clc 
% clear 
% close all; 
% 
% disp('Demo shows PSI method')
% disp('Simulated data habituate between stimulation with 0.1% factor')
% disp('Habituation continuous between trials'); 
% 
% NumStimulation = 100;
% NumTrials = 1; 
% grain     = 50;  % Out pdf is grain x grain = pdf = alpha x beta
% PM.PF = @LogisticFunc;
% StimulationResolution = 50; 
% 
% %parameter to simulate observer
% paramsGen = [10, 1, .02, .02]; 
% 
% %Stimulus values the method can select from
% PM.stimRange = (linspace(PM.PF([paramsGen(1) paramsGen(2) 0 0],.01,'inverse'),PM.PF([paramsGen(1) paramsGen(2) 0 0],.99,'inverse'),StimulationResolution));
% 
% %Define parameter ranges to be included in posterior
% priorAlphaRange = linspace(PM.PF([paramsGen(1) paramsGen(2) 0 0],.01,'inverse'),PM.PF([paramsGen(1) paramsGen(2) 0 0],.99,'inverse'),grain);
% priorBetaRange =  linspace(log10(.0625),log10(5),grain); %OBS. Stated in Log!
% % Jenny numbers
% %priorBetaRange =  linspace(log10(11.48),log10(1.995*10^114),grain); %OBS. Stated in Log!
% priorGammaRange = .02;  
% priorLambdaRange = .02; 
% 
% [PM.priorAlphas, PM.priorBetas, PM.priorGammas, PM.priorLambdas] = ndgrid(priorAlphaRange,priorBetaRange,priorGammaRange,priorLambdaRange);
% 
% %PDF
%     % "First, a prior probability distribution p0(lambda) for the 
%     % psychometric functions must be set up" [Kontsevich]
%     prior = ones(length(priorAlphaRange),length(priorBetaRange),length(priorGammaRange),length(priorLambdaRange));
%     prior = prior./numel(prior); 
%     %PM.pdf = prior;  
%   
% %LOOK UP TABEL (LUT)
%     % "Second, to speed up the method, a look-up table of conditional
%     % probabilities p(r|lambda,x) should be computed" [Kontsevich]
%     for a = 1:length(priorAlphaRange)
%         for b = 1:length(priorBetaRange) %OBS. Not calculated in log!
%             for g = 1:length(priorGammaRange)
%                 for L = 1:length(priorLambdaRange) 
%                     for sLevel = 1:length(PM.stimRange)
%                         %sandsyndligheden for at for korrekt response ved en given parameter sammensætning og intensitet.  
%                         PM.LUT(a,b,g,L,sLevel) = PM.PF([priorAlphaRange(a), 10.^priorBetaRange(b), priorGammaRange(g), priorLambdaRange(L)], PM.stimRange(sLevel));
%                     end
%                 end
%             end
%         end 
%     end
%     
%     clear a b g L sLevel 
%     clear grain StimulationResolution  
%     
%     disp('OBS! with contour plot, the trial nr is reset to zero '); 
%     doPlot = input('Do not plot (0), plot threshold (1), plot threshold PF and contour(2) ?: ');
% 
%     
% %% (PSI METHODE 2/3)
% 
%     clear newIntensityIndexPosition 
%     
%     if (doPlot) 
%         figure(1) 
%         xlim([1 NumTrials*NumStimulation])
%         ylim([min(PM.stimRange) max(PM.stimRange)])
%         for i = 1:NumTrials-1
%             figure(1)
%             hold on
%             plot([NumStimulation*i NumStimulation*i],[min(PM.stimRange) max(PM.stimRange)],':b')
%         end 
%         %set(gcf, 'Position',  [40, 3010, 1000, 600])
%         xlabel('Trial number') 
%         ylabel('Stimulus intensity')
%         legend('a','b')
%         
%         if(doPlot == 2)
%             NumTrials = 1; 
% 
%             figure(2) 
%             [X,Y] = meshgrid(priorAlphaRange, priorBetaRange);
%             %set(gcf, 'Position',  [1100, 580, 800, 400])
% 
%             figure(3) 
%             [X,Y] = meshgrid(priorAlphaRange, 10.^priorBetaRange);
%             StimLevelsFine = [min(PM.stimRange):(max(PM.stimRange) - min(PM.stimRange))./1000:max(PM.stimRange)];
%             Fit_correct = PM.PF(paramsGen, StimLevelsFine); 
%             %set(gcf, 'Position',  [1100, 50, 800, 400])
%             input('start?');
%         end
%     end 
%                   
%                 
% %% Simulate data and update method (PSI METHODE 3/3) 
% 
% 
% for CurrentTrialNum = 1:NumTrials
%     for CurrentTrialStim = 1:NumStimulation
%         
%     PM.x = []; 
%     PM.threshold = []; 
%     PM.pdf = prior; 
%     [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
%     [~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
%     PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
%     PM.x(1) = PM.xCurrent;
%     move = NumStimulation*(CurrentTrialNum-1);          
% 
%     while length(PM.x) <= NumStimulation
%         response = rand(1) < PM.PF([paramsGen(1)*1.001^(length(PM.x)+move) paramsGen(2) paramsGen(3) paramsGen(4)], PM.xCurrent);    %simulate observer
% 
%         %update PM based on response
%         PM = UpdateFunc(PM, response); 
%         % Test --> PMtest = PAL_AMPM_updatePM(PMtest,response);
%         responses(length(PM.x)-1) = response; 
% 
%         % plot? 
%         if (doPlot) 
%             figure(1) 
%             hold on;        
%             plot( 1+move:length(PM.x)-1+move, PM.threshold, 'b')
%             plot(length(PM.x)-1+move, paramsGen(1)*1.001^(length(PM.x)+move), '.', 'color','#B1B1B1', 'linewidth',0.1)
%             
%             figure(1) 
%             hold on; 
%             if response 
%                 plot(length(PM.x)+move,PM.x(end),'ok','MarkerFaceColor','k');
%             else 
%                 plot(length(PM.x)+move,PM.x(end),'ok');
%             end
% 
%             if(doPlot == 2)
%                 figure(2) 
%                 contour(X,Y, flip(flip(PM.pdf),2))        
%                 xlabel('threshold')
%                 ylabel('slope')
%                 title('PDF')
%                 grid on; 
% 
%                 figure(3) 
%                 Fit = PM.PF([PM.threshold(end), 10.^PM.slope(end), 0.02, 0.02], StimLevelsFine);
%                 clf
%                 hold on; 
%                 plot(StimLevelsFine,Fit_correct,'k-','linewidth',2);
%                 plot(StimLevelsFine,Fit,'g-','linewidth',2);
%                 title('Psychophysical Function')
% 
%                 for i = 2:length(responses)
%                     if responses(i) 
%                         dotSize = sum(PM.x(i) == PM.x(i:-1:1) & responses(i) == responses(i:-1:1));                                          
%                         plot(PM.x(i),responses(i),'ok','MarkerFaceColor','k', 'linewidth',dotSize);
%                     else
%                         dotSize = sum(PM.x(i) == PM.x(i:-1:1) & responses(i) == responses(i:-1:1));                                          
%                         plot(PM.x(i),responses(i),'ok','linewidth',dotSize);
%                     end 
%                 end 
%                 drawnow
%             end
%         end
%     end
%     end 
%     fprintf('Finish trial %2.0f. Threshold estimate %4.2f. Threshold at start and final state %4.2f : %4.2f \n', CurrentTrialNum, PM.threshold(end),paramsGen(1)*1.001^(1+move),paramsGen(1)*1.001^(length(PM.x)+move))
% end 
% 
% % running = movmean(PM.threshold,30);
% % 
% % figure;
% % plot(running)
% % hold on
% % plot(PM.threshold)

%% Setup (PSI METHODE 1/3) 
clc 
clear 
close all; 

disp('Demo shows PSI method')
disp('Simulated data habituate between stimulation with 0.1% factor')
disp('Habituation continuous between trials'); 

NumStimulation = 100;
NumTrials = 1; 
grain     = 50; 
PM.PF = @LogisticFunc;
StimulationResolution = 50; 

%parameter to simulate observer
paramsGen = [10, 1, .02, .02]; 

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
    clear grain StimulationResolution  
    
    disp('OBS! with contour plot, the trial nr is reset to zero '); 
    doPlot = input('Do not plot (0), plot threshold (1), plot threshold PF and contour(2) ?: ');

    
%% (PSI METHODE 2/3)

    clear newIntensityIndexPosition 
    
    if (doPlot) 
        figure(1) 
        xlim([1 NumTrials*NumStimulation])
        ylim([min(PM.stimRange) max(PM.stimRange)])
        for i = 1:NumTrials-1
            figure(1)
            hold on
            plot([NumStimulation*i NumStimulation*i],[min(PM.stimRange) max(PM.stimRange)],':b')
        end 
        %set(gcf, 'Position',  [40, 3010, 1000, 600])
        xlabel('Trial number') 
        ylabel('Stimulus intensity')
        
        if(doPlot == 2)
            NumTrials = 1; 

            figure(2) 
            [X,Y] = meshgrid(priorAlphaRange, priorBetaRange);
            %set(gcf, 'Position',  [1100, 580, 800, 400])

            figure(3) 
            [X,Y] = meshgrid(priorAlphaRange, 10.^priorBetaRange);
            StimLevelsFine = [min(PM.stimRange):(max(PM.stimRange) - min(PM.stimRange))./1000:max(PM.stimRange)];
            Fit_correct = PM.PF(paramsGen, StimLevelsFine); 
            %set(gcf, 'Position',  [1100, 50, 800, 400]) % Where the run
            %figure is placed on the computer
            input('start?');
        end
    end 
                  
                
%% Simulate data and update method (PSI METHODE 3/3) 

for CurrentTrialNum = 1:NumTrials
    % All this under is reseted during each trial
    PM.x = []; % Empty array to the stimuli intensities 
    PM.threshold = []; % Empty array to threshold values
    PM.pdf = prior; % Resets prior values 
    
    [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);
    % Caluclates the new intensity position 
    [~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);
    PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
    PM.x(1) = PM.xCurrent; % States that on the first array number
    move = NumStimulation*(CurrentTrialNum-1);          

    while length(PM.x) <= NumStimulation
        response = rand(1) < PM.PF([paramsGen(1)*1.001^(length(PM.x)+move) paramsGen(2) paramsGen(3) paramsGen(4)], PM.xCurrent);    %simulate observer

        %update PM based on response
        PM = UpdateFunc(PM, response); 
        % Test --> PMtest = PAL_AMPM_updatePM(PMtest,response);
        responses(length(PM.x)-1) = response; 

        % plot? 
        if (doPlot) 
            figure(1) 
            hold on;        
            plot( 1+move:length(PM.x)-1+move, PM.threshold, 'b','linewidth',2)
            plot(length(PM.x)-1+move, paramsGen(1)*1.001^(length(PM.x)+move), '.', 'color','#B1B1B1', 'linewidth',0.4)
            
            figure(1) 
            hold on; 
            if response 
                plot(length(PM.x)+move,PM.x(end),'ok','MarkerFaceColor','k');
            else 
                plot(length(PM.x)+move,PM.x(end),'ok');
            end

            if(doPlot == 2)
                % Contour plot
                figure(2) 
                contour(X,Y, flip(flip(PM.pdf),2))        
                xlabel('threshold')
                ylabel('slope')
                title('PDF')
                grid on; 

                % PF function plot
                figure(3) 
                Fit = PM.PF([PM.threshold(end), 10.^PM.slope(end), 0.02, 0.02], StimLevelsFine);
                clf
                hold on; 
                plot(StimLevelsFine,Fit_correct,'k-','linewidth',2);
                plot(StimLevelsFine,Fit,'b-','linewidth',2);
                title('Psychophysical Function')

                % Response and threshold fit-plot
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
    end
    fprintf('Finish trial %2.0f. Threshold estimate %4.2f. Threshold at start and final state %4.2f : %4.2f \n', CurrentTrialNum, PM.threshold(end),paramsGen(1)*1.001^(1+move),paramsGen(1)*1.001^(length(PM.x)+move))
end 

