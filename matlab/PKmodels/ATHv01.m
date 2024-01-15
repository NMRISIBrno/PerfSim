function curves = ATHv01(config, aif, t, parameters)
% t [min], aif [min]

    TRF = pkmodel([parameters.Fp, ...
                   parameters.E, ...
                   parameters.ve, ...
                   parameters.Tc], ...
                   t);

    %% Extract parameters from config and tissue structure
    NN = config.acquisition.timeAxis.NN;
    N = config.acquisition.timeAxis.N;
    R10 = 1 / parameters.T10;
    R20 = 1 / parameters.T2star0;
    
    %% Get concentration curves              
    conc = real(ifft(  fft(TRF,NN) .* fft(aif,NN)  )); % multiplying with Ts [min] done in pkmodel
    conc = conc(1:int32(N));    

    %% Conversion to relaxation rates
    R1 = conc * config.acquisition.r1 + R10; % R1(t) curve
    R2 = conc * parameters.r2star + R20;     % R2star(t) curve

    %% Create output table
    curves = table({TRF},{conc},{R1},{R2},...
        'VariableNames',{'TRF','ct','R1','R2'});
    curves.Properties.VariableUnits = {'1/min','mmol/l','1/s','1/s'};
    
end


function [TRF] = pkmodel(p,t)
% ATH model - equations from paper ...
%   Input:  p - vector of parameters [Fp E ve Tc]
%           t - time axis [min]
%   Output: TRF - tissue residue function

    Ts = t(2)-t(1); % sampling period       [min] 
    F = p(1);       % plasma flow           [ml/min/ml]
    E = p(2);       % extraction fraction   [-]
    ve = p(3);      % EES volume            [ml/ml]
    Tc = p(4);      % mean of transit times [min]
    sigma = Ts*2; %0.01; 
    Kep = E*F/ve;

    N = normcdf(0,Tc,sigma);
    erf2 = @(x)erf3(x,N);

    Hv = 1-(erf2((t-Tc)/(sqrt(2)*sigma))+erf2(Tc/(sqrt(2)*sigma)));

    Hp = E*exp(1/2*Kep^2*sigma^2+Kep*(Tc-t)).*(erf2((t-Tc)/(sqrt(2)*sigma)-Kep*sigma/sqrt(2))+erf2(Tc/(sigma*sqrt(2))+Kep*sigma/sqrt(2)));
    Hp(isnan(Hp)) = 0;
    h = Hp+Hv;
    h = h*F*Ts;
    TRF = h;
end

function h = erf3(x,N)
    y = sqrt(2)*x;
    h = (normcdf(y)-0.5)/(1-N);
end