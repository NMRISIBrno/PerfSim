function curves = TwoCX(config, aif, t, parameters)
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
% 2CXM model - equations from paper by Koh et al. 2011
%   Input:  p - vector of parameters [Fp E ve Tc]
%           t - time axis [min]
%   Output: TRF - tissue residue function

    Ts = t(2)-t(1);     % sampling period       [min] 
    Fp = p(1);          % plasma flow           [ml/min/ml]
    E = p(2);           % extraction fraction   [-]
    ve = p(3);          % EES volume            [ml/ml]
    Tc = p(4);          % mean of transit times [min]

    PS = Fp*E/(1-E);    % permeability-surface product [ml/min/ml] (for compartment model!)
    vp = Fp*Tc;         % plasma volume         [ml/ml]
    
    alpha = 0.5*(-(PS/vp + PS/ve + Fp/vp) + sqrt((PS/vp + PS/ve + Fp/vp).^2 - 4*(PS*Fp)/(ve*vp))); 
    beta = 0.5*(-(PS/vp + PS/ve + Fp/vp) - sqrt((PS/vp + PS/ve + Fp/vp).^2 - 4*(PS*Fp)/(ve*vp))); 
    A = (alpha + PS/vp + PS/ve) ./ (alpha-beta);
    
    h = A * exp(alpha.*t) + (1-A) * exp(beta.*t); 
    h = h*Fp*Ts;
    TRF = h;
    
end
