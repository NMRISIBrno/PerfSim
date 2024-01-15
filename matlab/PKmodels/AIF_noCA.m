function [aif] = AIF_noCA(p, Ts, N, delay)

% Input:  p     - structure with AIF parameters [min]!
%         Ts    - sampling period [s]! (conversion to mins here in fcn)
%         N     - number of samples
%         delay - bolus arrival time [s]! 
% Output: aif_triexponencial - modeled aif with delay [mM = mmol/l]


%% Create AIF for precontrast scans (vector of zeros)

t = (0:(N-1))*Ts/60; % time axis in minutes
aif = zeros(size(t));


end