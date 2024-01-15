function [SI, SITE0] = IRLL(config, R1, R2, TEind)
% R1, R2 in [1/s]

    krho = config.acquisition.krho;
    FA = config.acquisition.FA; 
    
    tau = config.acquisition.tau;   % time between excitations         
    td = config.acquisition.td;     % time between IR and first excitation
    tr = config.acquisition.tr;     % relaxation time at the end of excitation train
    nProjPerInv = config.acquisition.nProjPerInv;

    Etau = exp(-tau.*R1);
    Etd = exp(-td.*R1);
    Etr = exp(-tr.*R1);
    
    F = (1-Etau)./(1-cos(FA)*Etau);
    Q = (-F.*cos(FA).*Etr.*Etd.*(1-(cos(FA).*Etau).^(nProjPerInv-1))-2.*Etd+Etr.*Etd+1)./(1+cos(FA).*Etr.*Etd.*(cos(FA).*Etau).^(nProjPerInv-1));
    
    SI = zeros(1,config.acquisition.kSampling.projections);
    i = 1;
    while i <= config.acquisition.kSampling.projections
        for n=1:config.acquisition.nProjPerInv
            SI(i) = sin(FA)*(krho+0*1i)*(F(i)+(cos(FA).*Etau(i)).^(n-1).*(Q(i)-F(i)));
            i = i+1;
        end
    end

    SITE0 = [];

end
