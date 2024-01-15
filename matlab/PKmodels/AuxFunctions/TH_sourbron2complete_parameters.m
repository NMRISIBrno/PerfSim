function [param_complete,K,label_out]=TH_sourbron2complete_parameters(par_matrix)
% [parametry,K,label_out]=TH_sourbron2complete_parameters(par_matrix)
% input:
%   par_matrix=[Fp T Tc Te tau] - original values
% output:
%   param_complete=[F_blood F_plasma T Tc Te E PS v_p v_e K_trans Kep tau]
%   K - transformation matrix for extension of original covariance
%       matrix, i.e. COVB_complete=K*COVB*K'
%   label_out={'F_b' 'F_p' 'T' 'T_c' 'T_e' 'E' 'PS' 'v_p' 'v_e' 'K^t^r^a^n^s' 'k_e_p' '\tau'}

Hct=0.42;
param_complete=zeros(size(par_matrix,1),12);

if nargout>1
    K=zeros([12,5,size(par_matrix,1)]);
end

if nargout>2
    label_out={'F_b' 'F_p' 'T' 'T_c' 'T_e' 'E' 'PS' 'v_p' 'v_e' 'K^t^r^a^n^s' 'k_e_p' '\tau'};
end

for n=1:size(par_matrix,1)
    par=par_matrix(n,:);
    Fp=par(1);
    T=par(2);
    Tc=par(3);
    Te=par(4);
    tau=par(5);
    
    F_blood=Fp/(1-Hct);
    E=1-exp(-(T-Tc)/Te);
    PS=Fp*(T-Tc)/Te;    
    vp=Fp*Tc;
    ve=Fp*(T-Tc);
    Ktrans=Fp*E;
    kep=Ktrans/ve;
    
    param_complete(n,:)=[F_blood Fp T Tc Te E PS vp ve Ktrans kep tau];
    
    if nargout>1
        % par_matrix=[Fp T Tc Te tau] - original values
        % par_out=[F_blood F_plasma T Tc Te E PS v_p v_e K_trans Kep tau]
        K(1,1)=1/(1-Hct); %[dFb/dFp]
        K(2,1)=1; %[dFp/dFp]
        K(3,2)=1; %[dT/dT]
        K(4,3)=1; %[dTc/dTc]
        K(5,4)=1; %[dTe dTe]
        K(6,2:4)=[exp(-(T-Tc)/Te)/Te -exp(-(T-Tc)/Te)/Te -exp(-(T-Tc)/Te)*(T-Tc)/Te^2]; %[dE/d(T Tc Te)]
        K(7,1:4)=[(T-Tc)/Te, Fp/Te -Fp/Te, -Fp*(T-Tc)/Te^2];% [dPS/d(Fp T Tc Te)]
        K(8,[1 3])=[Tc Fp];%[dvp/d(Fp Tc)]
        K(9,1:3)=[T-Tc Fp -Fp];% [dve/d(Fp T Tc)]
        K(10,1:4)=[E, Fp/Te*exp(-(T-Tc)/Te), -Fp/Te*exp(-(T-Tc)/Te), -Fp*(T-Tc)/Te^2*exp(-(T-Tc)/Te)];%[dKtrans/d(Fp T Tc Te)]
        K(11,2)=-E/(T-Tc)^2+exp(-(T-Tc)/Te)/(Te*(T-Tc));% [dkep/dT]
        K(11,3:4)=[-K(11,2), -exp(-(T-Tc)/Te)/Te^2];% [dkep/d(Tc Te)]
        K(12,5)=1; %[dtau/dtau]
    end
end