function config = GenerateAIF(config)
% GENERATEAIF Generate AIF based on model given in config.
%   Output file: 'AIF.mat' - vector of AIF in time

    fh = str2func(config.phantom.AIF.model);
    
    %% Extract sampling parameters
    Ts = config.acquisition.timeAxis.Ts; % sampling period [s]
    N = config.acquisition.timeAxis.N;   % number of samples
    
    %% Generate AIF curve (delayed)
    aif  = fh(config.phantom.AIF,... 
              Ts,... 
              N,...
              config.acquisition.delay);
    
    %% Save AIF and update config
    config.phantom.AIF.fileName = 'AIF.mat';
    save([config.workDir.output, config.phantom.AIF.fileName], 'aif');

end