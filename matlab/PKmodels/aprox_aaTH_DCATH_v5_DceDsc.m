function curves = aprox_aaTH_DCATH_v5_DceDsc(config, aif, t, parameters)
% t [min], aif [min]
    
    [TRF,Hp,He] = pkmodel([parameters.Fp, ...
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
    cp = real(ifft(  fft(Hp,NN) .* fft(aif,NN)  )); 
    cp = cp(1:int32(N));
    ce = real(ifft(  fft(He,NN) .* fft(aif,NN)  )); 
    ce = ce(1:int32(N));

    %% GCM model for DceDsc
    fh = str2func(parameters.GCM_model{1,1});
    R2 = fh([parameters.r2p, parameters.r2e],cp,ce) + R20; % R2star(t) curve
    
    %% Conversion to relaxation rates
    R1 = conc * config.acquisition.r1 + R10; % R1(t) curve

    %% Create output table
    curves = table({TRF},{conc},{cp},{ce},{R1},{R2},...
        'VariableNames',{'TRF','ct','cp','ce','R1','R2'});
    curves.Properties.VariableUnits = {'1/min','mmol/l','mmol/l','mmol/l','1/s','1/s'};
    
end

function [h,Hp,He] = pkmodel(x,t)
% DCATH model (aaJW model with dispersed Tc), paper by Koh et al.
%   -The inclusion of capillary distribution in the adiabatic tissue 
%    homogeneity model of blood flow
% h = DCATH_v5(x,t), where x=[F E ve Tc sigma]

    Ts = t(2)-t(1); % sampling period       [min]
    F = x(1);       % plasma flow           [ml/min/ml]
    E = x(2);       % extraction fraction   [-]
    ve = x(3);      % EES volume            [ml/ml]
    Tc = x(4);      % mean of transit times [min] (not sensitive - is multiplied by 10 to increase influence)  
    sigma = Ts*2;   % std.deviation of transit times    

    Kep = E*F/ve;

    N = normcdf(0,Tc,sigma);
    erf2 = @(x)erf3(x,N);
    
    Hp = 1-(erf2((t-Tc)/(sqrt(2)*sigma))+erf2(Tc/(sqrt(2)*sigma)));

    He = E*exp(1/2*Kep^2*sigma^2+Kep*(Tc-t)).*(erf2((t-Tc)/(sqrt(2)*sigma)-Kep*sigma/sqrt(2))+erf2(Tc/(sigma*sqrt(2))+Kep*sigma/sqrt(2)));
    He(isnan(He)) = 0;
    
    h = He+Hp;
    h = h*F*Ts;
    
    % Model for plasma (p) and (ees - e) compartments concentrations
    Hp = (Hp*F*Ts)/(F*Tc); % vp = (F*Tc)
    He = (He*F*Ts)/ve;

end

function h = erf3(x,N)
    y = sqrt(2)*x;
    h = (normcdf(y)-0.5)/(1-N);
end