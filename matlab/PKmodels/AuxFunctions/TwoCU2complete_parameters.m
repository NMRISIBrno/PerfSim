function [param_complete,covariance_matrices,label_out]=TwoCU2complete_parameters(par_matrix,covariance_matrices,x_fixed)
% [parametry,K,label_out]=TwoCU2complete_parameters(par_matrix)
% input:
%   par_matrix=[Fp Tp E tau] - original values or
%   par_matrix=[Fb Fp Tp Tc E PS vp Ktrans tau] - complete parameters
% output:
%   param_complete=[F_blood F_plasma Tp Tc E PS v_p K_trans tau]
%   param_complete=[Fp Tp E tau] - parameters as input to the related file,
%   i.e.: TwoCU_FT.m
%   K - transformation matrix for extension of original covariance
%       matrix, i.e. COVB_complete=K*COVB*K'
%   label_out={'F_b' 'F_p' 'T_p' 'T_c' 'E' 'PS' 'v_p' 'K^t^r^a^n^s' '\tau'}

switch size(par_matrix,2)
    case 4
        Hct=0.42;
        param_complete=zeros(size(par_matrix,1),9);
        if nargout>1
            covariance_matrices=zeros([9,4,size(par_matrix,1)]); % Hack
        end
        if (nargout>1&&nargin>1)
            %             K=zeros([12,5,size(par_matrix,1)]);
            K=zeros(9,4);
            if ~iscell(covariance_matrices)
                covariance_matrices={covariance_matrices};
            end
            if (size(par_matrix,2)~=size(covariance_matrices{1},2))
                covariance_matrices=cellfun(@(x)covmtx_fixed(x,x_fixed),covariance_matrices,'UniformOutput',false);
            end
        end
        
        if nargout>2
            label_out={'F_b' 'F_p' 'T_p' 'T_c' 'E' 'PS' 'v_p' 'K^t^r^a^n^s' '\tau'};
        end
        
        for n=1:size(par_matrix,1)
            par=par_matrix(n,:);
            Fp=par(1);
            Tp=par(2);
            E=par(3);
            tau=par(4);
            
            F_blood=Fp/(1-Hct);
            Ktrans=Fp*E;
            PS=Ktrans/(1-E);
            vp=Tp*(Fp+PS);
            Tc=vp/Fp;
            
            param_complete(n,:)=[F_blood Fp Tp Tc E PS vp Ktrans tau];
            
            if (nargout>1&&nargin>1)
                % par_matrix=[Fp Tp E tau] - original values
                % par_out=[F_blood F_plasma Tc E PS v_p K_trans tau]
                K(1,1)=1/(1-Hct); %[dFb/dFp]
                K(2,1)=1; %[dFp/dFp]
                K(3,2)=1; %[dTp/dTp]
                K(4,2:3)=[1 Tp/(1-E)]/(1-E);%[dTc/d(Tp E)]
                K(5,3)=1; %[dE/dE]
                K(6,[1 3])=[E Fp-Fp*E/(1-E)]/(1-E); %[dPS/d(Fp E)]
                K(7,1:3)=[Tp Fp Tp*Fp/(1-E)]/(1-E);%[dvp/d(Fp Tp E)]
                K(8,[1 3])=[E Fp];% [dK_trans/d(Fp E)]
                K(9,4)=1; %[dtau/dtau]
                
                covariance_matrices{n}=K*covariance_matrices{n}*K.';
            end
        end
    case 9
        label_out={'F_p' 'T_p' 'E' '\tau'};
        param_complete=par_matrix(:,[2 3 5 9]);
        
        if (nargout>1&&nargin>1) % tohle tady k nicemu neni
            K=zeros(9,4);
            K(2,1)=1;
            K(3,2)=1;
            K(5,3)=1;
            K(9,5)=1;
        end
    otherwise
        error('Wrong number of the parameters in the input matrix.')
end
end
function [K]=covmtx_fixed(K,x_fixed)
I=zeros(length(x_fixed));
I(~x_fixed,~x_fixed)=K;
K=I;
end