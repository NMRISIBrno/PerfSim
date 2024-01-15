function curves = TwoCU_FT(config, aif, t, parameters)


    %% Extract parameters from config and tissue structure            
    R10 = 1 / parameters.T10;
    R20 = 1 / parameters.T2star0;
    if verLessThan('matlab','9.4')
        Ts = parameters.Properties.UserData.SamplingPeriod; % sampling period [s]
    else
        Ts = parameters.Properties.CustomProperties.SamplingPeriod; % sampling period [s]
    end
    
    %% PK model rendered in Fourier domain
    % aif - Sampled
    global aif_database
    aif_database = struct('model',aif,'t_cut',t(end),'Ts',Ts,'x_aif',[],'index',zeros(1,1,'uint16'),'data',[]);
    model.aif = struct('model',@aif_nonblind_real,'no_of_par',1);
    x_aif_label = 'fake';
    no_of_aif_par = 1;            

    model.irf = struct('model',@pkmodelTwoCU_FT,'no_of_par',5);
%                     x_irf_label={'F_p' 'T_p' 'E' '\tau'};
    x_label = {'fake' 'F_p' 'T_p' 'E' '\tau' 't_trunc'};
    rescale_parameters = ones(1,model.irf.no_of_par+1);
    x_fixed = true(1,model.irf.no_of_par+1);
%   Recompute perametrisation:  [Fp E ve Tc] -> x_irf_label based on Sourbron [2013]
    x = [0 parameters.Fp parameters.Tc*(1-parameters.E) parameters.E 0 t(end)];
    x(x==0) = eps;

%     % rof GT only:
%     %x_model_roi(roi,:) = x;
%     % get the complete set of perfusion parameters in x_full (used previously as GT)
%     [x_full, K, x_irf_label_all] = TwoCU2complete_parameters(x(:,(no_of_aif_par+1):end-1));
%     x_full = [x(:,1:no_of_aif_par) x_full];
    
    [conc, dCfit, aif_and_h] = C_fit_FT(x,t,model,rescale_parameters,x_fixed);
    TRF = aif_and_h.irf;
    aif = aif_and_h.aif;
    conc = conc.';            

    %% Conversion to relaxation rates
    R1 = conc * config.acquisition.r1 + R10; % R1(t) curve
    R2 = conc * parameters.r2star + R20;     % R2star(t) curve

    %% Create output table
    curves = table({aif},{TRF},{conc},{R1},{R2},...
        'VariableNames',{'AIF','TRF','ct','R1','R2'});
    curves.Properties.VariableUnits = {'-','-','mmol/l','l/s','l/s'};
    
end


function [H,dH] = pkmodelTwoCU_FT(par,w,deriv)
% parametrization [Fp Tp E tau]
% 2CU model, from Sourbron 2012 Table A.2,
% derivatives from TwoCU_symbolic.m

w=w(:);
Fp=par(1);
T_p=par(2);
E=par(3);
tau=par(4);
t_cut=par(5);

phase_shift=exp(-1j*par(4)*w);

p=w*1j;
% s=p;

% part to be truncated, i.e. shifted GK (Kety) model
H_GK=1./(p+1/T_p);
H_GK=H_GK.*exp(-1/T_p.*(t_cut-tau))*(1-E); % reduce energy
% H_GK=H_GK.*exp(-p*(t_cut-tau)); % shift to t_cut

% % part to be truncated, flat area
H_flat=E./p; % GK model
H_flat(1)=tau*E;
% H_flat=H_flat.*exp(-p*(t_cut-tau)); % shift to t_cut

% Truncation together
H_GK = H_GK + H_flat;
H_GK=H_GK.*exp(-p*(t_cut-tau)); % shift to t_cut

% 2CU model, from Sourbron 2012 Table A.2:
H_2CU = (1-E)*T_p./(1+p*T_p) + E./p;
H_2CU(1) = (1-E)*T_p + E*t_cut; % mean value set manually

H=H_2CU-H_GK; % truncation
H=H.*phase_shift; % time-shift

if nargout>1 % if derivatives needed
    if nargin>2
        ind=1:length(deriv);
        ind=ind(deriv);
        dH=zeros(length(w),nnz(deriv));
        dH_GK=zeros(length(w),nnz(deriv));
        for n=1:length(ind)
            switch ind(n)
                case 1
                    % dH/dFp
                    dH(:,n)=H./Fp;
                case 2
                    % dH/dTp
                    dH(:,n)=-(E-1.0)./(T_p.*p+1.0)+T_p.*p.*1.0./(T_p.*p+1.0).^2.*(E-1.0);
                    dH_GK(:,n)=-(exp(-(t_cut-tau)./T_p).*1.0./(T_p.*p+1.0).^2.*(E-1.0).*(T_p+t_cut-tau+T_p.*p.*t_cut-T_p.*p.*tau))./T_p;
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 3
                    % dH/dE
                    % dH is defined with rect pulse, dH_GK is just
                    % exponential
                    dH(:,n)=-T_p./(T_p.*p+1.0)+1./p.*(1-exp((tau-t_cut).*p));
                    dH(1,n)=-T_p-(tau-t_cut); % for p=0 it is minus sum of Fp and mean of rect pulse, i.e. 1*pulse_width
                    dH_GK(:,n)=-(T_p.*exp(-(t_cut-tau)./T_p))./(T_p.*p+1.0).*exp(-1j*(t_cut-tau)*w);
                    dH(:,n)=dH(:,n)-dH_GK(:,n);
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 4
                    % dH/dtau
                    dH(:,n)=H.*-1j.*w; % time derivative, neni zkonrolovano
                    dH_GK(:,n)=exp(-((T_p.*p+1.0).*(t_cut-tau))./T_p)-E.*exp(-((T_p.*p+1.0).*(t_cut-tau))./T_p)+E.*exp(-p.*(t_cut-tau));
                    dH(:,n)=dH(:,n)+dH_GK(:,n).*phase_shift;
            end
        end
    else
        dH=zeros(length(w),5);
        dH(:,1)=H./Fp; % together with dH_GK/dF
        dH(:,2)=-(E-1.0)./(T_p.*p+1.0)+T_p.*p.*1.0./(T_p.*p+1.0).^2.*(E-1.0);
        dH(:,3)=[-T_p-(tau-t_cut), -T_p./(T_p.*p(2:end)+1.0)+1./p(2:end).*(1-exp((tau-t_cut).*p(2:end)))];
        
        dH(:,1)=0; % it is done different way
        dH_GK(:,2)=-(exp(-(t_cut-tau)./T_p).*1.0./(T_p.*p+1.0).^2.*(E-1.0).*(T_p+t_cut-tau+T_p.*p.*t_cut-T_p.*p.*tau))./T_p;
        dH_GK(:,3)=-(T_p.*exp(-(t_cut-tau)./T_p))./(T_p.*p+1.0).*exp(-1j*(t_cut-tau)*w);
        
        dH(:,4)=H.*-1j.*w; % time derivative, neni zkonrolovano
        dH_GK(:,4)=exp(-((T_p.*p+1.0).*(t_cut-tau))./T_p)-E.*exp(-((T_p.*p+1.0).*(t_cut-tau))./T_p)+E.*exp(-p.*(t_cut-tau));
        dH(:,4)=dH(:,4)+dH_GK(:,4).*phase_shift;
        
        dH=dH(:,2:3)-dH_GK(:,2:3).*repmat(exp(-1j*(t_cut-tau)*w),1,3); % truncation of derivatives
        
        dH=dH(:,2:3).*repmat(phase_shift,1,2); % time-shift of derivatives
        
    end
    dH=Fp.*dH;
end
H=Fp.*H;



end