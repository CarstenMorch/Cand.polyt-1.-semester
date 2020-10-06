%
%LogisticFunction   Estimate values on Logistic Psychometric Function
%
%Syntax: 
%   y = LogisticFunc(parameter, x)
%
%Inputs: 
%   'parameter': vector containing the Psychometric Funtion coefficients 
%   (i.e., [alpha, beta, gamma, lambda]). If gamma and lambda is not 
%   specificed, they will by default be set to zero.    
%
%   'x': vector containing values of which the values of the Psychometric 
%   funtion shall be evaluated. 
%
%Output:
%   y: values of Psychometric function evaluated at values in 'x'
%
%Example: 
%   StimLevels = [0:0.5:10]; 
%   pcorrect = LogisticFunction([3 1 0.1 0], StimLevels); 
%
%latest edited: 
%   06-10-2020

function y = LogisticFunc(parameter,x)

gamma = 0;
lambda = 0;
alpha = parameter(1);
beta = parameter(2);
if length(parameter) > 2
    gamma = parameter(3);
    if length(parameter) > 3
        lambda = parameter(4);
    end
end

y = gamma + (1 - gamma - lambda).*(1./(1+exp(-1*(beta).*(x-alpha))));



