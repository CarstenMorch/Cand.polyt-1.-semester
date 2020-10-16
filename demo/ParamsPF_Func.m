%% Function to unpack the parameters
% With inspiration in PAL_unpackParamsPF from the Palemedes toolbox
% Input parameters are 'params' which is a struct
% The function assigns values to alpha, beta, gamma and lambda 
function [alpha, beta, gamma, lambda] = ParamsPF_Func(params)
    % gamma and lambda are set to a fixed value if not set to other earlier
    gamma = 0; % This should be 0.02
    lambda = 0; % This should be 0.02
    
    if isstruct(params) %returns 1 if the struct 'params' is true, if not, continues to 'else'
        alpha = params.alpha; % assigns alpha from params to value alpha
        beta = params.beta; % assigns beta from params to value beta
        if isfield(params,'gamma')  % returns false if any of the fields are empty 
            gamma = params.gamma; %assigns gamma the value of params.gamma
            if isfield(params,'lambda') % checks that lambda is a field
                lambda = params.lambda; % assigns params.lambda value to lambda
            end 
        end
    else % If 'isstruct(params)' returns false, don't exsist, it is made here
        alpha = params(1); % the first spot in the 'params' struct is assigned to alpha
        beta = params(2); % The second 'params' spot is assigned to beta 
        if length(params) > 2 % If the struct is longer than 2 
            gamma = params(3); % Gamma is assigned to 'params' struct 3rd spot
            if length(params) > 3 % If the struct is longer than 3
                lambda = params(4); % Lambda is assigned to the 4th spot in 'params'
            end
        end 
    end
