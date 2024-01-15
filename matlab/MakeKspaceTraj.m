function config = MakeKspaceTraj(config)
% MAKEKSPACETRAJ Generate kspace trajectory with chosen sampling and save
% time axis parameters into config.
%   Output file: 'kTrajectory.mat' - matrices kx and ky with trajectory


    %% Get parameters from config and generate trajectory
    nSamples = config.acquisition.kSampling.nSamples;   % number of samples of each echo (sampling rate is fixed, given by the FOV)

    % Create acquisition scheme:
    % list of acquisition steps, for each TR (i.e. for each echo, or set of
    % echoes for array coils) the "acquired" echo signal is specified as one
    % entry in the following two vectors:
    %
    %   angles - angle of the acquisition trajectory line [rad]
    %   phEncShifts - shift in the phase encoding direction [samples in the k-space]

    switch config.acquisition.kSampling.method
        case 'cartesian'
            % Cartesian sampling (for Michal Bartos, or for reference):
            phEncSteps  = config.acquisition.kSampling.phEncSteps; 
            repetitions = config.acquisition.kSampling.repetitions; 
            angles = zeros(1, phEncSteps*repetitions);      % angles - for Cartesian sampling always zero
            phEncShifts = flip((-phEncSteps/2):(phEncSteps/2-1)); % phase encoding step increments
            phEncShifts = repmat(phEncShifts, 1, repetitions);
            
            % k-space coordinates
            kx = ( flip((-nSamples/2):(nSamples/2-1)) )' * cos(angles);
            ky = ( flip((-nSamples/2):(nSamples/2-1)) )' * sin(angles) + repmat(phEncShifts,[nSamples, 1]);
            
            % normalization to [-0.5, 0.5]
            kx = kx/nSamples;
            ky = ky/nSamples;
            
            % second normalization to 0.5*[-Nsamples/xsize, Nsamples/ysize] 
            kx = kx*nSamples/config.phantom.xSize;
            ky = ky*nSamples/config.phantom.ySize;
            
        case 'radial'
            % Radial 2D sampling with set angle increment:
            projections = config.acquisition.kSampling.projections; 
            angleIncrement = config.acquisition.kSampling.angleIncrement;
            angles = (0:(projections-1)) * angleIncrement; 
            phEncShifts = zeros(1,projections); % phase encoding shifts - for radial sampling always zero
            
            % k-space coordinates
            kx = ( (-nSamples/2):(nSamples/2-1) )' * cos(angles);
            ky = ( (-nSamples/2):(nSamples/2-1) )' * sin(angles) + repmat(phEncShifts,[nSamples, 1]);
            
            % normalization to [-0.5, 0.5]
            kx = kx/nSamples;
            ky = ky/nSamples;
            
            % second normalization to 0.5*[-Nsamples/xsize, Nsamples/ysize]
            kx = kx*nSamples/config.phantom.xSize;
            ky = ky*nSamples/config.phantom.ySize;

        case 'rosettes'
            rosettes = config.acquisition.kSampling.rosettes; 
            angleIncrement = config.acquisition.kSampling.angleIncrement;
            angles = (0:(rosettes-1)) * angleIncrement; 
            w1 = config.acquisition.kSampling.w1; %[rad/s]
            w2 = config.acquisition.kSampling.w2; %[rad/s]
            deltaT = 1/config.acquisition.rBW; % [s]
            t = (0:(nSamples-1)) * deltaT;
            t = t(:);
            kmax = 0.5;
            k = kmax * ( sin(w1*t) .* exp(1i*w2*t) ) * exp(1i*angles); % nSamples x numberOfAngles (numberOfExcitations)
            kx = real(k);
            ky = imag(k);
            
            % normalization from [-0.5, 0.5] to 0.5*[-config.acquisition.kSampling.Nres/xsize, config.acquisition.kSampling.Nres/ysize]
            kx = kx*config.acquisition.kSampling.Nres/config.phantom.xSize;
            ky = ky*config.acquisition.kSampling.Nres/config.phantom.ySize;

    end
    
    %% Save k-space trajectory and update config
    config.acquisition.kSampling.fileName = 'kTrajectory.mat';
    save([config.workDir.output, config.acquisition.kSampling.fileName], 'kx', 'ky');
    
end