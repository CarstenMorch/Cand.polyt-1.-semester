% Bias calculation
% Based on Kontsevich bias of the threshold estimate after j-th trial

function y = bias_calc(thresholdRand,thresholdEst,NumExpCur,NumExpTotal,NumStimCur,NumStimTotal)
    for i = NumExpCur : NumExpTotal
        for j = NumStimCur:NumStimTotal
           s =  sum(log(thresholdRand)-log(thresholdEst));
           y = (s/NumExpTotal) * 20*log10(20); % should be dB
        end
    end
end