%% The psychometric function
% With inspiration from PAL_Gumbel from the Palmedes toolbox
% Input: params and varargin(which collects all inputs)
function y = PsychmetricFunc(params,x,varargin)
% Calling the 'ParamsPF_Func' funtion, which makes 'params' a struct
[alpha, beta, gamma, lambda] = ParamsPF_Func(params);

if ~isempty(varargin) % runs if varargin is not an empty array
    if strncmpi(varargin{1},'Inverse',3)% Compares the first 3 elements of varargin{1} to the inverse
        c = (x-gamma)./(1-gamma-lambda)-1;
        c = -1.*log(-1.*c);
        c = log10(c);
        y = alpha + c./beta;
    end
    if strncmpi(varargin{1},'Derivative',3)
        y = (1-gamma-lambda).*exp(-1.*10.^(beta.*(x-alpha))).*log(10).*10.^(beta.*(x-alpha)).*beta;
    end 
else
    y = gamma+(1-gamma-lambda).*(1-exp(-1.*10.^(beta.*(x-alpha))));
end
