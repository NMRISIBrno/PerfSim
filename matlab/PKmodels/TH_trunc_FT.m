function curves = TH_trunc_FT(config, aif, t, parameters)
% t - time in [min]

%% Extract parameters from config and tissue structure
R10 = 1 / parameters.T10;
R20 = 1 / parameters.T2star0;
Ts = parameters.Properties.CustomProperties.SamplingPeriod; % sampling period [s]

%% PK model rendered in Fourier domain
% aif - Sampled
global aif_database
aif_database = struct('model',aif,'t_cut',t(end),'Ts',Ts,'x_aif',[],'index',zeros(1,1,'uint16'),'data',[]);
model.aif = struct('model',@aif_nonblind_real,'no_of_par',1);
x_aif_label = 'fake';
no_of_aif_par = 1;

model.irf = struct('model',@pkmodelTH_trunc_FT,'no_of_par',6);
rescale_parameters = ones(1,model.irf.no_of_par+1);
x_fixed = true(1,model.irf.no_of_par+1);
x_label = {'fake' 'Fp' 'vp' 've' 'PS' 'tau' 't_trunc'};
%                   Recompute parametrisation:  [Fp E ve Tc] -> x_irf_label
x = [0 parameters.Fp, parameters.Fp*parameters.Tc, parameters.ve ...
    -parameters.Fp*log(1-max(parameters.E,eps)), 0, t(end)];
x(x==0) = eps;

%     % rof GT only:
%     x_model_roi(roi,:) = x;
%     % get the complete set of perfusion parameters in x_full (used previously as GT)
%     [x_full, K, x_irf_label_all] = TH_garpebring2complete_parameters(x(:,(no_of_aif_par+1):end-1));
%     x_full = [x(:,1:no_of_aif_par), x_full];

[conc, dCfit, aif_and_h] = C_fit_FT(x,t,model,rescale_parameters,x_fixed);
TRF = aif_and_h.irf;
aif = aif_and_h.aif;
conc = conc.';

%% Conversion to relaxation rates
R1 = conc * config.acquisition.r1 + R10; % R1(t) curve
R2 = conc * parameters.r2star + R20;     % R2star(t) curve

%% Create output table
curves = table({TRF},{conc},{R1},{R2},...
    'VariableNames',{'TRF','ct','R1','R2'});
curves.Properties.VariableUnits = {'1/min','mmol/l','1/s','1/s'};

end



function [H,dH] = pkmodelTH_trunc_FT(par, w, deriv)
% [H dH]=TH_FT([Fp vp ve PS tau],w,(deriv))
% computes Tissue Homogeneity (TH) model in frequency domain based on Garpebring 2009 paper:
% optional deriv has the form deriv=[true false true true true]
w=w(:);
Fp=par(1);
vp=par(2);
ve=par(3);
PS=par(4);
tau=par(5);
t_cut=par(6);
% tau=exp(-1j*par(5)*(0:(length(w)-1)));
phase_shift=exp(-1j*par(5)*w);

a=PS/Fp;
Tc=vp/Fp;
c=ve/vp;
E=(1-exp(-a));
Ktrans=Fp*E;

p=w*1j;

% part to be truncated, i.e. shifted GK (Kety) model
H_GK=1./(p/Ktrans+1/ve); % GK model
H_GK=H_GK*exp(-Ktrans/ve*(t_cut-tau-Tc)); % reduce energy
H_GK=H_GK.*exp(-1j*(t_cut-tau)*w); % shift to t_cut

H_TH=(exp(-(a+Tc*p))-1).*(a+Tc*p).*(a*ve+vp*(p*c*Tc+a))./...
    (a^2*(exp(-(a+Tc*p))-1)-(a+Tc*p).*p*Tc.*(c*(a+Tc*p)+a));
% H=H*sqrt(length(w));

H=H_TH-H_GK; % truncation
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
                    dH(:,n)=-((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)-((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).*(vp.*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*ve)+1.0./Fp.^2.*PS.*ve))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)+(exp(-PS./Fp-(p.*vp)./Fp).*(PS./Fp+(p.*vp)./Fp).*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)-1.0./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp).^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp).*(1.0./Fp.^3.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*-2.0+1.0./Fp.^2.*PS.^2.*exp(-PS./Fp-(p.*vp)./Fp).*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp)+1.0./Fp.^2.*p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp)+(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp))./Fp+(p.*vp.*(1.0./Fp.^2.*PS+(ve.*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp);
                    %                     dH_GK(:,n)=-exp(-p.*t_cut).*exp((Fp.*(exp(-PS./Fp)-1.0).*(t_cut-vp./Fp))./ve).*1.0./(1.0./ve-p./(Fp.*(exp(-PS./Fp)-1.0))).^2.*((1.0./Fp.^2.*p)./(exp(-PS./Fp)-1.0)+1.0./Fp.^3.*PS.*p.*exp(-PS./Fp).*1.0./(exp(-PS./Fp)-1.0).^2)+(exp(-p.*t_cut).*exp((Fp.*(exp(-PS./Fp)-1.0).*(t_cut-vp./Fp))./ve).*(((exp(-PS./Fp)-1.0).*(t_cut-vp./Fp))./ve+(vp.*(exp(-PS./Fp)-1.0))./(Fp.*ve)+(PS.*exp(-PS./Fp).*(t_cut-vp./Fp))./(Fp.*ve)))./(1.0./ve-p./(Fp.*(exp(-PS./Fp)-1.0)));
                    %                     dH_GK(:,n)=(-(1.0./Fp.^2.*Ktrans.*vp.*exp((Ktrans.*(-t_cut+par(5)+vp./Fp))./ve))./(ve.*(p./Ktrans+1.0./ve))).*(1-exp(-a));
                    %                     dH_GK(:,n)=((exp((Ktrans.*(-t_cut+tau+vp./Fp))./ve).*1.0./(Ktrans+p.*ve).^2.*(Ktrans.^2.*vp-Fp.*Ktrans.^2.*t_cut+Fp.*Ktrans.^2.*tau+Fp.*p.*ve.^2+Ktrans.*p.*ve.*vp-Fp.*Ktrans.*p.*t_cut.*ve+Fp.*Ktrans.*p.*tau.*ve))./Fp).*(-(PS.*exp(-PS./Fp))./Fp+1.0);
                    dH_GK(:,n)=-(exp(-(Fp.*(exp(-PS./Fp)-1.0).*(-t_cut+tau+vp./Fp))./ve).*(((exp(-PS./Fp)-1.0).*(-t_cut+tau+vp./Fp))./ve-(vp.*(exp(-PS./Fp)-1.0))./(Fp.*ve)+(PS.*exp(-PS./Fp).*(-t_cut+tau+vp./Fp))./(Fp.*ve)))./(1.0./ve-p./(Fp.*(exp(-PS./Fp)-1.0)))-exp(-(Fp.*(exp(-PS./Fp)-1.0).*(-t_cut+tau+vp./Fp))./ve).*1.0./(1.0./ve-p./(Fp.*(exp(-PS./Fp)-1.0))).^2.*((1.0./Fp.^2.*p)./(exp(-PS./Fp)-1.0)+1.0./Fp.^3.*PS.*p.*exp(-PS./Fp).*1.0./(exp(-PS./Fp)-1.0).^2);
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 2
                    dH(:,n)=((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*ve)./Fp).*(PS./Fp+(p.*vp)./Fp))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)+1.0./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp).^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp).*((p.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp+1.0./Fp.^3.*PS.^2.*p.*exp(-PS./Fp-(p.*vp)./Fp)+1.0./Fp.^2.*p.^2.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp)-(p.*vp.*(ve.*1.0./vp.^2.*(PS./Fp+(p.*vp)./Fp)-(p.*ve)./(Fp.*vp)).*(PS./Fp+(p.*vp)./Fp))./Fp)+(p.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(Fp.*(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp))-(p.*exp(-PS./Fp-(p.*vp)./Fp).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(Fp.*(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp));
                    % %                     dH_GK(:,n)=-(exp(-p.*t_cut).*exp((Fp.*(exp(-PS./Fp)-1.0).*(t_cut-vp./Fp))./ve).*(exp(-PS./Fp)-1.0))./(ve.*(1.0./ve-p./(Fp.*(exp(-PS./Fp)-1.0))));
                    %                     dH_GK(:,n)=(Ktrans.*exp((Ktrans.*(-t_cut+par(5)+vp./Fp))./ve))./(Fp.*ve.*(p./Ktrans+1.0./ve));
                    dH_GK(:,n)=(Ktrans.*exp((Ktrans.*(-t_cut+tau+vp./Fp))./ve))./(Fp.*ve.*(p./Ktrans+1.0./ve));
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 3
                    dH(:,n)=((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).^2)./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)+(p.*1.0./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp).^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).^3.*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./Fp;
                    %                     dH_GK(:,n)=(Fp.*exp(-p.*t_cut).*exp(-PS./Fp).*exp((Fp.*(exp(-PS./Fp)-1.0).*(t_cut-vp./Fp))./ve).*(exp(PS./Fp)-1.0).^2.*1.0./(-Fp+Fp.*exp(PS./Fp)+p.*ve.*exp(PS./Fp)).^2.*(Fp.*vp-Fp.^2.*t_cut+Fp.^2.*t_cut.*exp(PS./Fp)+Fp.*ve.*exp(PS./Fp)-Fp.*vp.*exp(PS./Fp)-p.*ve.*vp.*exp(PS./Fp)+Fp.*p.*t_cut.*ve.*exp(PS./Fp)))./ve;
                    %                     dH_GK(:,n)=(Fp.*exp(-p.*t_cut).*exp(-PS./Fp).*exp((Fp.*(exp(-PS./Fp)-1.0).*(t_cut-vp./Fp))./ve).*(exp(PS./Fp)-1.0).^2.*1.0./(-Fp+Fp.*exp(PS./Fp)+p.*ve.*exp(PS./Fp)).^2.*(Fp.*vp-Fp.^2.*t_cut+Fp.^2.*t_cut.*exp(PS./Fp)+Fp.*ve.*exp(PS./Fp)-Fp.*vp.*exp(PS./Fp)-p.*ve.*vp.*exp(PS./Fp)+Fp.*p.*t_cut.*ve.*exp(PS./Fp)))./ve;
                    %                     dH_GK(:,n)=-(Ktrans.^2.*exp((Ktrans.*(-t_cut+par(5)+vp./Fp))./ve).*1.0./(Ktrans+p.*ve).^2.*(-Fp.*ve+Ktrans.*vp-Fp.*Ktrans.*t_cut+Fp.*Ktrans.*par(5)+p.*ve.*vp-Fp.*p.*t_cut.*ve+Fp.*p.*par(5).*ve))./(Fp.*ve);
                    dH_GK(:,n)=-(Ktrans.^2.*exp((Ktrans.*(-t_cut+tau+vp./Fp))./ve).*1.0./(Ktrans+p.*ve).^2.*(-Fp.*ve+Ktrans.*vp-Fp.*Ktrans.*t_cut+Fp.*Ktrans.*tau+p.*ve.*vp-Fp.*p.*t_cut.*ve+Fp.*p.*tau.*ve))./(Fp.*ve);
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 4
                    dH(:,n)=((ve./Fp+vp./Fp).*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)+((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(Fp.*(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp))+1.0./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp).^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp).*(1.0./Fp.^2.*PS.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*-2.0+1.0./Fp.^3.*PS.^2.*exp(-PS./Fp-(p.*vp)./Fp)+1.0./Fp.^2.*p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp)+(p.*vp.*(1.0./Fp+ve./(Fp.*vp)).*(PS./Fp+(p.*vp)./Fp))./Fp)-(exp(-PS./Fp-(p.*vp)./Fp).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(Fp.*(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp));
                    %                     dH_GK(:,n)=(Fp.*exp(-p.*t_cut).*exp(-PS./Fp).*exp((Fp.*(exp(-PS./Fp)-1.0).*(t_cut-vp./Fp))./ve).*(exp(PS./Fp)-1.0).^2.*1.0./(-Fp+Fp.*exp(PS./Fp)+p.*ve.*exp(PS./Fp)).^2.*(Fp.*vp-Fp.^2.*t_cut+Fp.^2.*t_cut.*exp(PS./Fp)+Fp.*ve.*exp(PS./Fp)-Fp.*vp.*exp(PS./Fp)-p.*ve.*vp.*exp(PS./Fp)+Fp.*p.*t_cut.*ve.*exp(PS./Fp)))./ve;
                    %                     dH_GK(:,n)=((exp((Ktrans.*(-t_cut+par(5)+vp./Fp))./ve).*1.0./(Ktrans+p.*ve).^2.*(Ktrans.^2.*vp-Fp.*Ktrans.^2.*t_cut+Fp.*Ktrans.^2.*par(5)+Fp.*p.*ve.^2+Ktrans.*p.*ve.*vp-Fp.*Ktrans.*p.*t_cut.*ve+Fp.*Ktrans.*p.*par(5).*ve))./Fp).*exp(-a);
                    %                     dH_GK(:,n)=((exp((Ktrans.*(-t_cut+tau+vp./Fp))./ve).*1.0./(Ktrans+p.*ve).^2.*(Ktrans.^2.*vp-Fp.*Ktrans.^2.*t_cut+Fp.*Ktrans.^2.*tau+Fp.*p.*ve.^2+Ktrans.*p.*ve.*vp-Fp.*Ktrans.*p.*t_cut.*ve+Fp.*Ktrans.*p.*tau.*ve))./Fp).*exp(-a);
                    dH_GK(:,n)=(exp(-PS./Fp).*exp(-(Fp.*(exp(-PS./Fp)-1.0).*(-t_cut+tau+vp./Fp))./ve).*(-t_cut+tau+vp./Fp))./(ve.*(1.0./ve-p./(Fp.*(exp(-PS./Fp)-1.0))))+1.0./Fp.^2.*p.*exp(-PS./Fp).*exp(-(Fp.*(exp(-PS./Fp)-1.0).*(-t_cut+tau+vp./Fp))./ve).*1.0./(exp(-PS./Fp)-1.0).^2.*1.0./(1.0./ve-p./(Fp.*(exp(-PS./Fp)-1.0))).^2;
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 5
                    dH(:,n)=H_TH*-1j.*w; % time derivative
                    %                     dH_GK(:,n)=Ktrans*exp(-Ktrans/ve*(t_cut-tau-Tc)).*exp(-1j*(t_cut-tau-Tc)*w);
                    %                     dH_GK(:,n)=Ktrans*exp(-Ktrans/ve*(t_cut-tau-Tc)).*exp(-1j*(t_cut-tau)*w);
                    dH_GK(:,n)=-(Fp.*exp(-(Fp.*(exp(-PS./Fp)-1.0).*(-t_cut+tau+vp./Fp))./ve).*(exp(-PS./Fp)-1.0))./(ve.*(1.0./ve-p./(Fp.*(exp(-PS./Fp)-1.0))));
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                    
            end
        end
    else
        dH=zeros(length(w),5);
        dH(:,1)=-((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)-((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).*(vp.*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*ve)+1.0./Fp.^2.*PS.*ve))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)+(exp(-PS./Fp-(p.*vp)./Fp).*(PS./Fp+(p.*vp)./Fp).*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)-1.0./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp).^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp).*(1.0./Fp.^3.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*-2.0+1.0./Fp.^2.*PS.^2.*exp(-PS./Fp-(p.*vp)./Fp).*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp)+1.0./Fp.^2.*p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp)+(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp))./Fp+(p.*vp.*(1.0./Fp.^2.*PS+(ve.*(1.0./Fp.^2.*PS+1.0./Fp.^2.*p.*vp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp);
        dH(:,2)=((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*ve)./Fp).*(PS./Fp+(p.*vp)./Fp))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)+1.0./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp).^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp).*((p.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp+1.0./Fp.^3.*PS.^2.*p.*exp(-PS./Fp-(p.*vp)./Fp)+1.0./Fp.^2.*p.^2.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp)-(p.*vp.*(ve.*1.0./vp.^2.*(PS./Fp+(p.*vp)./Fp)-(p.*ve)./(Fp.*vp)).*(PS./Fp+(p.*vp)./Fp))./Fp)+(p.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(Fp.*(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp))-(p.*exp(-PS./Fp-(p.*vp)./Fp).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(Fp.*(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp));
        dH(:,3)=((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).^2)./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)+(p.*1.0./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp).^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).^3.*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./Fp;
        dH(:,4)=((ve./Fp+vp./Fp).*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp))./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp)+((exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(Fp.*(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp))+1.0./(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp).^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp).*(1.0./Fp.^2.*PS.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0).*-2.0+1.0./Fp.^3.*PS.^2.*exp(-PS./Fp-(p.*vp)./Fp)+1.0./Fp.^2.*p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp)+(p.*vp.*(1.0./Fp+ve./(Fp.*vp)).*(PS./Fp+(p.*vp)./Fp))./Fp)-(exp(-PS./Fp-(p.*vp)./Fp).*(PS./Fp+(p.*vp)./Fp).*(vp.*(PS./Fp+(p.*ve)./Fp)+(PS.*ve)./Fp))./(Fp.*(1.0./Fp.^2.*PS.^2.*(exp(-PS./Fp-(p.*vp)./Fp)-1.0)-(p.*vp.*(PS./Fp+(ve.*(PS./Fp+(p.*vp)./Fp))./vp).*(PS./Fp+(p.*vp)./Fp))./Fp));
        
        dH_GK(:,1)=(p/Ktrans+1/ve).^-2 .*p/Ktrans^2 * -PS/(Fp*exp(a));
        dH_GK(:,3)=(p/Ktrans+1/ve).^-2 *ve^-2;
        dH_GK(:,4)=(p/Ktrans+1/ve).^-2 .*p/Ktrans^2 * exp(-a);
        
        dH=dH-dH_GK; % truncation of derivatives
        
        dH=dH(:,1:4).*repmat(phase_shift,1,4); % time-shift of derivatives
        dH(:,5)=H_TH*-1j.*w; % time derivative
    end
end


end