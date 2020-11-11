
%% Load data 
clc
clear 
close all; 

load('PM1201_09.mat')
load('resp1201_09.mat')

%% Brute force search (Maximum likelihood 1/2)

PF = @LogisticFunc;
stimIntst = PM.x(1:end-1); 
%NumPos = [59 53 68 83 92 99]; %number of trials with correct response.
%OutOfNum = [100 100 100 100 100 100]; %number of stimulus in each stimLevel

%Define sg for brute force search 
sg.alpha = linspace(4,7,100);
sg.beta = 1; %linspace(log10(.0625),log10(5),100); 
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
                temp(stimNum) = log(PF([sg.alpha(a)+sg.delta(d)*stimNum, sg.beta, sg.gamma, sg.lambda], stimIntst(stimNum)));
            else
                temp(stimNum) = log(1-PF([sg.alpha(a)+sg.delta(d)*stimNum, sg.beta, sg.gamma, sg.lambda], stimIntst(stimNum)));
            end 
        end
        likelihoodGrid(a,d)= sum(temp); 
    end
end 
[M,I] = max(likelihoodGrid(:)); 
[I_row, I_col] = ind2sub(size(likelihoodGrid),I)
parameters = [sg.alpha(I_row) sg.delta(I_col)]



% contour(likelihoodGrid)
% 
% % Find maximum value in likelihood grid 
% Amax = max(likelihoodGrid(:));                                           % Maximum Value
% Idx = find(likelihoodGrid(:) == Amax);                                   % Returns Linear Indices To All ‘Amax’ Values
% [r,c,dim3,dim4] = ind2sub(size(likelihoodGrid), Idx);
% % Return found values 
% parameters = [sg.alpha(r) sg.beta(c) sg.gamma(dim3) sg.lambda(dim4)];


% likelihoodGrid = zeros(length(sg.alpha),length(sg.beta),length(sg.gamma),length(sg.lambda),length(sg.delta));    
% for a = 1:length(sg.alpha)
%     for b = 1:length(sg.beta)
%         for d = 1:length(sg.delta)
%             for g = 1:length(sg.gamma)
%                 for L = 1:length(sg.lambda) 
%                     for stimNum = 1:length(stimIntst)
%                         temp = zeros(1,length(stimIntst));  
%                         if responses(stimNum) == true
%                             temp(stimNum) = PF([sg.alpha(a)+sg.delta*stimNum, sg.beta(b),sg.gamma(g), sg.lambda(L)], stimIntst(stimNum));
%                         else
%                             temp(stimNum) = 1-PF([sg.alpha(a)+sg.delta*stimNum, sg.beta(b), sg.gamma(g), sg.lambda(L)], stimIntst(stimNum));
%                         end 
%                         likelihoodGrid(a,b,g,L,stimNum)=prod(temp); 
%                         clear temp; 
%                     end
%                 end
%             end
%         end
%     end 
% end 
