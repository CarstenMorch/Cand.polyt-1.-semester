%
%LogisticFunction   Estimate Logistic Psychometric Function values
%
%   syntax: y = LogisticFunction(parameter, x)
%
%   y = LogisticFunction(parameter,x), where 'parameter' contains the four
%   parameters of a Psychometric Funtion (i.e., [alpha beta gamma lambda]),
%   returns the Psychometric Function evaluated at values in 'x'. 'x' is
%   array of any size.
%
% Created by Buus, 30-09-2020
% version 1 

function y = LogisticFunction(parameter,x)

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



