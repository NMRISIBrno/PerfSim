function aif = AIF_triexpG(p, Ts, N, delay)
% AIF_TRIEXPG Triexponential AIF model. 
%   Parameters described in DOI:10.1109/EMBC.2014.6944569.
%   Input:  p     - structure with AIF parameters [min]!
%           Ts    - sampling period [s]! (conversion to mins here in fcn)
%           N     - number of samples
%           delay - bolus arrival time [s]! 
%   Output: aif - modeled aif with delay [mM = mmol/l]


%% Sum of 2 decreasing exponentials and one gamma-variate function
t = (0:(N-1))*Ts/60; % time axis in minutes

s = zeros(2,int32(N));

s(1,:) = t.^p.beta .* p.A.*(exp(-p.tau1*t));
s(2,:) = t.^p.beta .* p.B.*(exp(-p.tau2*t));
s(3,:) = t.^p.beta .* p.C.*(exp(-p.tau3*t));

aif = sum(s);

delaySamples = round(delay / Ts);
aif  = [zeros(1,delaySamples), aif(1:(N-delaySamples))]; % delay AIF
   
end