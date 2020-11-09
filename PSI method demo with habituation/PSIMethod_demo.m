
%% Setup (PSI METHODE 1/3) 
clc 
clear 
close all; 

disp('Demo shows PSI method')
disp('Simulated data habituate between stimulation with 0.1% factor')
disp('Habituation continuous between trials'); 

% Sets up parameters
NumStimulation = 40; % Number of stimulation in each trial
NumTrials = 6; % Number of trials
grain     = 50; % Resolution of the parameters during the code
StimulationResolution = 50; % Resolution of the stimulation

% The psykometric function - saves in a struct "PM"
PM.PF = @LogisticFunc;

%parameter to simulate observer
paramsGen = [10, 1.05, .02, .02]; 

%Stimulus values the method can select from
PM.stimRange = (linspace(PM.PF([paramsGen(1) paramsGen(2) 0 0],.01,'inverse'),PM.PF([paramsGen(1) paramsGen(2) 0 0],.99,'inverse'),StimulationResolution));

%Define parameter ranges to be included in posterior
priorAlphaRange = linspace(PM.PF([paramsGen(1) paramsGen(2) 0 0],.01,'inverse'),PM.PF([paramsGen(1) paramsGen(2) 0 0],.99,'inverse'),grain);
priorBetaRange =  linspace(log10(.0625),log10(5),grain); %OBS. Stated in Log!
priorGammaRange = .02;  % Fixed value
priorLambdaRange = .02; % Fixed value

[PM.priorAlphas, PM.priorBetas, PM.priorGammas, PM.priorLambdas] = ndgrid(priorAlphaRange,priorBetaRange,priorGammaRange,priorLambdaRange);

%PDF
    % "First, a prior probability distribution p0(lambda) for the 
    % psychometric functions must be set up" [Kontsevich]
    prior = ones(length(priorAlphaRange),length(priorBetaRange),length(priorGammaRange),length(priorLambdaRange)); % Construcs matrix of 0, of the size of the parameters
    prior = prior./numel(prior); % Total number of elements in the prior
    PM.pdf = prior;   % Saves the prior in the PM struct 
  
%LOOK UP TABEL (LUT)
    % "Second, to speed up the method, a look-up table of conditional
    % probabilities p(r|lambda,x) should be computed" [Kontsevich]
   % Develops a 'for'-loop that loops through the length of the parameters
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
    
%% Deletes unessary variables 
    clear a b g L sLevel 
    clear StimulationResolution newIntensityIndexPosition 
    
    disp('OBS. Trial is set to 1, if 2 is entered');
    doPlot = input('Do not plot (0), plot threshold (1), plot threshold PF and contour(2) ?: ');

    
%% (PSI METHODE 2/3)

% Creates an 'if' loop, that plots the data
    if (doPlot) 
        figure(1) % Fit to the threshold and responses  
        xlim([1 NumTrials*NumStimulation])
        ylim([min(PM.stimRange) max(PM.stimRange)])
        
        for i = 1:NumTrials-1
            figure(1)
            hold on
            plot([NumStimulation*i NumStimulation*i],[min(PM.stimRange) max(PM.stimRange)],':b')
        end 
        %set(gcf, 'Position',  [40, 3010, 1000, 600])
        xlabel('Stimulation number') 
        ylabel('Stimulus intensity')
        
        if(doPlot == 2)
            NumTrials = 1; 
            figure(2) % pdf contour plot
            [X,Y] = meshgrid(priorAlphaRange, priorBetaRange);
            %set(gcf, 'Position',  [1100, 580, 800, 400])

            figure(3) % Psychometric function
            [X,Y] = meshgrid(priorAlphaRange, 10.^priorBetaRange);
            StimLevelsFine = [min(PM.stimRange):(max(PM.stimRange) - min(PM.stimRange))./1000:max(PM.stimRange)];
            Fit_correct = PM.PF(paramsGen, StimLevelsFine); 
            %set(gcf, 'Position',  [1100, 50, 800, 400])
            input('Press enter to start');
        end
    end 
                  
                
%% Simulate data and update method (PSI METHODE 3/3) 

for CurrentTrialNum = 1:NumTrials 
    PM.x = [];  % Empty array to the stimuli intensities 
    PM.threshold = [];  % Empty array to threshold values
    PM.pdf = prior;  % Resets prior values 
    
% Step 1 to 2
    % Uses the function 'PosteriorNextTrialFunc', and puts the data ind the PM
    % struct - Uses the pdf and LUT to calculate the outcomes
    [PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure,PM.pSuccessGivenx] = PosteriorNextTrailFunc(PM.pdf, PM.LUT);

% Step 3 to 5 
    % Uses the function 'EntropyFunc', to calculate the new intensity index
    % position. ~ is the logical 'not'
    [~, newIntensityIndexPosition] = EntropyFunc(PM.PosteriorNextTrailSuccess,PM.PosteriorNextTrialFailure, PM.pSuccessGivenx);

% Next stimulation intensity - calculates the stim range, from the new
% intensity
    PM.xCurrent = PM.stimRange(newIntensityIndexPosition);
    PM.x(1) = PM.xCurrent; 
    
% Moving window length (MARIA: not quite sure of this)
% Needs an explantion of this one...  
move = NumStimulation*(CurrentTrialNum-1); 

% Creates a while loop, that loops until the length is not longer than the
% number for trials
    while length(PM.x) <= NumStimulation
    % Simulates the observer response - returns 1 or 0, if the response is
    % greater (1) or smaller than the output from the PF struct
        response = rand(1) < PM.PF([paramsGen(1)*1.00001^(length(PM.x)+move) paramsGen(2) paramsGen(3) paramsGen(4)], PM.xCurrent);    %simulate observer
        fprintf('Trial %4.f. current stimulus level %4.2f. Response %1.0f \n', length(PM.x), PM.xCurrent , response)

        responses(length(PM.x)) = response; 

    % Updates the PM based on the response, this is done by the function
    % UpdateFunc
        PM = UpdateFunc(PM, response); 

       % Plots data
        if (doPlot) 
            figure(1) %Fit to the threshold and responses  
            hold on;        
            plot( 1+move:length(PM.x)-1+move, PM.threshold, 'b') % The fitted blue line
            plot(length(PM.x)-1+move, paramsGen(1)*1.00001^(length(PM.x)+move), '.', 'color','#B1B1B1', 'linewidth',0.1) % The theoretically assmued threshold
            
            % Plot of the responses
            figure(1) 
            hold on; 
            if response 
                plot(length(PM.x)-1+move,PM.x(end-1),'ok','MarkerFaceColor','k');
            else 
                plot(length(PM.x)-1+move,PM.x(end-1),'or');
            end

            if(doPlot == 2)
                figure(2) % pdf Contourplot
                contour(X,Y, flip(flip(PM.pdf),2))     % flippes order of the elements   
                xlabel('threshold')
                ylabel('slope')
                title('PDF')
                grid on; 

                figure(3) % PF 
                Fit = PM.PF([PM.threshold(end), 10.^PM.slope(end), 0.02, 0.02], StimLevelsFine);
                clf
                hold on; 
                plot(StimLevelsFine,Fit_correct,'k-','linewidth',2); % The theroretical PF
                plot(StimLevelsFine,Fit,'g-','linewidth',2); % The fitted PF
                title('Psychophysical Function')

                % Plots the responses on fig 3 - PF
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

