function config = AdjustAcqLen(config)

    % Adjust acquisition length
    truncScanTime = config.acquisition.kSampling.repetitions * (config.acquisition.timeAxis.Ts * config.acquisition.kSampling.phEncSteps);
    N = ceil(truncScanTime / config.acquisition.timeAxis.Ts); 
    Nmin = N * 2 - 1; 
    NN = (2^ceil(log2(Nmin))); 

    % Save into config
    config.acquisition.timeAxis.scanTime = truncScanTime;
    config.acquisition.timeAxis.N = N;
    config.acquisition.timeAxis.NN = NN;
    
end


