function [SI, SITE0] = FLASH2D(config, R1, R2, TEind)
% R1, R2 in [1/s]

    krho = config.acquisition.krho;
    FA = config.acquisition.FA;
    TR = config.acquisition.TR;
    TE = config.acquisition.TE;

    SITE0 = krho*sin(FA) * (1-exp(-R1*TR)) ./ (1-cos(FA)*exp(-R1*TR));
    SI = SITE0 .* exp(-TE(TEind)*R2);

end